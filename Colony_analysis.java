package LeidenUniv.ColonyMorphology;
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import LeidenUniv.Tools.*;
import ij.plugin.frame.*;
import java.util.*;
import ij.io.*;
import java.io.*;
import ij.measure.*;

/**
 * @author willemsejj
 * @version 1.1.0
 * 
 */
public class Colony_analysis implements PlugIn {
	FileSaver fs;
	int imageSize;
	public String imageTitle;
	/**
	 * @param arg no arguments are taken into this plugin
	 */
	public void run(String arg) {
		ImagePlus image = IJ.getImage();
		imageTitle=image.getTitle().substring(0,image.getTitle().lastIndexOf("."));
		String ImageDir = image.getOriginalFileInfo().directory; // ends with a \
		String outputdir = ImageDir+image.getTitle().substring(0,image.getTitle().lastIndexOf("."));
		new File(outputdir).mkdir();
		GenericDialog gd = null;
		gd = new GenericDialog("Analysis settings");
		gd.enableYesNoCancel("I already have an stack of images to analyze", "I want to mark colonies from a scanned plate");
		gd.addNumericField("Estimated colony size in pixels", 150,0);
		gd.showDialog();

		if (gd.wasCanceled()){
	            IJ.log("User clicked 'Cancel', no data is generated");
			} else if (gd.wasOKed()) {
				imageSize=(int)gd.getNextNumber();
				//image.hide();
				getData(image,outputdir);
				IJ.log("Colonies are analysed, data is stored in "+outputdir);
			} else {
				imageSize=(int)gd.getNextNumber();
				ImagePlus Icol = getColonyStack (image,outputdir);
				if (Icol!=null){
					getData(Icol,outputdir);
					IJ.log("Colonies are analysed, data is stored in "+outputdir);
				} else {
					IJ.log("No colonies selected, no data is generated");
				}
			}
		
	}
	
	/**
	 * @param image The original image from which the colonies are sselected
	 * @param outputdir The directory where data is saved
	 * @return The stack of images taken from getColonies
	 */
	public ImagePlus getColonyStack (ImagePlus image, String outputdir){
		// this gets the selected colonies
		// we need some alignment regions around the colonies
		// thus the area selected needs to be bigger than the actual colonies
		ImagePlus Icol = getColonies(image, imageSize);
		if (Icol!=null){IJ.log("Colonies are analysed, data is stored in "+outputdir);
			fs = new FileSaver(Icol);
			//outputdir+"/selected_colonies.tif"
			String cFile=outputdir+"/selected_colonies.tif";
			File cout = new File(cFile);
			int count=0;
			while (cout.exists()){
				count++;
				cFile=outputdir+"/selected_colonies_"+count+".tif";
				cout= new File(cFile);
			}
			fs.saveAsTiff(cFile);
			return Icol;
		} else {
			return null;
		}
	}
	
	/**
	 * @param Icol The imagestack containing all colonies to be analysed
	 * @param outputdir The directory where data is saved
	 */
	public void getData (ImagePlus Icol, String outputdir){
		
		//Icol.show();
		ImagePlus ori = Icol.duplicate();
		//ori.show();
		//ori.setTitle("this is original");
		
		

		// on the original image we are going to measure the shapes
		double [][] colData =getColonyInformation(ori, outputdir);	
		// coldata holds for each colony 29 parameters
		getColonyHistograms(ori, outputdir);
		// now we need to align them
		ImagePlus [] AI = getAlignedImages(Icol);
		Icol.hide();
		Icol.close(); //original image is no longer needed
		// now we have aligned images, 0 is the original data aligned, 1 is the edges aligned.
		// on 0 we need to do the radial transform and get pearson values for that
		// also we need to get the area, skew of all three colors, and roundness
		// on 1 we need to test the colocalisations
		ImagePlus RA = AI[0];
		fs = new FileSaver(RA);
		String cFile=outputdir+"/"+imageTitle+"aligned_colonies.tif";
		File cout = new File(cFile);
		int count=0;
		while (cout.exists()){
			count++;
			cFile=outputdir+"/"+imageTitle+"aligned_colonies_"+count+".tif";
			cout= new File(cFile);
		}
		fs.saveAsTiff(cFile);

		//AI[0].show();
		// this is an image stack, and for each layer in the stack we need to perform
		// radialtransfrom, averagefill, getcolorprofile
		// before we need to get the three colors per slice into different datasets.
		int [] dim = RA.getDimensions();
		int stdepth=Math.max(dim[3],dim[4]);
		ProfilePlot [][] profiledata = new ProfilePlot[stdepth][];
		double [][] Maxdata = new double [stdepth][];
		for (int i=0;i<stdepth;i++){
			ImagePlus rad = new ImagePlus(""+i,RA.getImageStack().getProcessor(i+1));
			ImagePlus imp3 = getRadialTransform(rad);
			imp3 = getHorizontalAverageFilled(imp3);
			//imp.show();
			ProfilePlot [] pfarray = getColorProfiles(imp3);
			profiledata [i]=pfarray;
			imp3.setTitle("why not");
			IJ.run(imp3, "RGB Color", "");
			IJ.run(imp3, "Gaussian Blur...", "sigma=3");
			Maxdata[i]=getMaxProfile(imp3);
		}

		// profiledata holds for each colony 3 profileplots (r,g,b)
		// Maxdata holds for each colony 1 double []

		// now we have the data, we want to calculate pearson correlation coefficients on all profiles
		// Furthermore we need to save all the data in a subdirectory of the image directory where the original image is located
		
		ImagePlus CA =AI[1];
		fs = new FileSaver(CA);
		cFile=outputdir+"/"+imageTitle+"aligned_edges.tif";
		cout = new File(cFile);
		count=0;
		while (cout.exists()){
			count++;
			cFile=outputdir+"/"+imageTitle+"aligned_edges_"+count+".tif";
			cout= new File(cFile);
		}
		fs.saveAsTiff(cFile);
		//CA.show();
		// this is an image stack, we need to analyze the colocalization on the stack. which can be done by splitting into 
		// three colors and analysing the seperate stacks

		// after we have all the data we need to create output files
		
		ImagePlus [] edgesplit = null;

		if (CA.getBitDepth()==16){
			edgesplit = new ImagePlus[3];
			Duplicator dc = new Duplicator();
			edgesplit[0]=dc.run(CA,1,1,1,dim[3],1,1);
			edgesplit[1]=dc.run(CA,2,2,1,dim[3],1,1);
			edgesplit[2]=dc.run(CA,3,3,1,dim[3],1,1);
			
		} else if (CA.getBitDepth()==24){
			edgesplit = new ChannelSplitter().split(CA);
		}	
		java.util.List<Double> dataRe = new ArrayList<Double>();
		java.util.List<Double> dataGe = new ArrayList<Double>();
		java.util.List<Double> dataBe = new ArrayList<Double>();
		java.util.List<Double> dataGre = new ArrayList<Double>();
		java.util.List<Double> dataRp = new ArrayList<Double>();
		java.util.List<Double> dataGp = new ArrayList<Double>();
		java.util.List<Double> dataBp = new ArrayList<Double>();
		java.util.List<Double> dataMp = new ArrayList<Double>();
		
		for (int i =0;i<dim[3]*dim[4]-1;i++){ // LOOP TRROUGH ALL COLONIES
			for (int j=i+1;j<dim[3]*dim[4];j++){ // LOOP ONLY THROUGH CONSECUTIVE COLONIES
				// edge comparisons
				float[][] i1 = edgesplit[0].getImageStack().getProcessor(i+1).getFloatArray();
				float[][] i2 = edgesplit[0].getImageStack().getProcessor(j+1).getFloatArray();
				float[][] i3 = edgesplit[1].getImageStack().getProcessor(i+1).getFloatArray();
				float[][] i4 = edgesplit[1].getImageStack().getProcessor(j+1).getFloatArray();
				float[][] i5 = edgesplit[2].getImageStack().getProcessor(i+1).getFloatArray();
				float[][] i6 = edgesplit[2].getImageStack().getProcessor(j+1).getFloatArray();
				float[][] i7 = CA.getImageStack().getProcessor(i+1).getFloatArray();
				float[][] i8 = CA.getImageStack().getProcessor(j+1).getFloatArray();
				// now we need to invoke the 2d getpearson
				dataRe.add(getPearson2d(i1,i2));
				dataGe.add(getPearson2d(i3,i4));
				dataBe.add(getPearson2d(i5,i6));			
				dataGre.add(getPearson2d(i7,i8));			

				// profile comparisons
				dataRp.add(getPearson1D(profiledata[i][0].getProfile(), profiledata[j][0].getProfile()));
				dataGp.add(getPearson1D(profiledata[i][1].getProfile(), profiledata[j][1].getProfile()));
				dataBp.add(getPearson1D(profiledata[i][2].getProfile(), profiledata[j][2].getProfile()));
				dataMp.add(getPearson1D(Maxdata[i], Maxdata[j]));
			}
		}
		// and now for writing the output
		// one file for the raw data
		// one file for all comparable data

		try {
			// write the raw data
			String rawFile=outputdir+"/"+imageTitle+"RawData.csv";
			File rawout = new File(rawFile);
			count=0;
			while (rawout.exists()){
				count++;
				rawFile=outputdir+"/"+imageTitle+"RawData_"+count+".csv";
				rawout= new File(rawFile);
			}
			PrintWriter out = new PrintWriter(rawFile);
			out.println("\"sep=,\"");
			out.println("Red profiles");	

			for (int i=0;i<profiledata.length;i++){
				out.print("Colony"+(i+1)+",");
				String pdata = Arrays.toString(profiledata[i][0].getProfile());
				out.println(pdata.substring(1,pdata.length()-1));		
			}

			out.println("");	
			out.println("Green profiles");	
			for (int i=0;i<profiledata.length;i++){
				out.print("Colony"+(i+1)+",");
				String pdata = Arrays.toString(profiledata[i][1].getProfile());
				out.println(pdata.substring(1,pdata.length()-1));		
			}
			
			out.println("");	
			out.println("Blue profiles");	
			for (int i=0;i<profiledata.length;i++){
				out.print("Colony"+(i+1)+",");
				String pdata = Arrays.toString(profiledata[i][1].getProfile());
				out.println(pdata.substring(1,pdata.length()-1));		
			}

			out.println("");	
			out.println("Max difference profiles");	
			for (int i=0;i<Maxdata.length;i++){
				out.print("Colony"+(i+1)+",");
				String pdata = Arrays.toString(Maxdata[i]);
				out.println(pdata.substring(1,pdata.length()-1));		
			}
			
			out.close();
		}catch (Exception e){
			IJ.log(e.toString());
			StackTraceElement[] trace = e.getStackTrace();
			for (int i=0;i<trace.length;i++){
				IJ.log(trace[i].toString());
			}
		}
		try {
			//write the pearson correlations
			String pearFile=outputdir+"/"+imageTitle+"PearsonData.csv";
			File pearout = new File(pearFile);
			count=0;
			while (pearout.exists()){
				count++;
				pearFile=outputdir+"/"+imageTitle+"PearsonData_"+count+".csv";
				pearout= new File(pearFile);
			}
			PrintWriter out2 = new PrintWriter(pearFile);
			out2.println("\"sep=,\"");
			out2.println("Red Edge Correlations");
			out2.println(getPearsonTable(dataRe,Maxdata.length));
			out2.println("Green Edge Correlations");
			out2.println(getPearsonTable(dataGe,Maxdata.length));
			out2.println("Blue Edge Correlations");
			out2.println(getPearsonTable(dataBe,Maxdata.length));
			out2.println("Grey Edge Correlations");
			out2.println(getPearsonTable(dataGre,Maxdata.length));
			out2.println("Red Profile Correlations");
			out2.println(getPearsonTable(dataRp,Maxdata.length));
			out2.println("Green Profile Correlations");
			out2.println(getPearsonTable(dataGp,Maxdata.length));
			out2.println("Blue Profile Correlations");
			out2.println(getPearsonTable(dataBp,Maxdata.length));
			out2.println("Maximum radial deviation Correlations");
			out2.println(getPearsonTable(dataMp,Maxdata.length));

			
			out2.close();
		}catch (Exception e){
			IJ.log(e.toString());
			StackTraceElement[] trace = e.getStackTrace();
			for (int i=0;i<trace.length;i++){
				IJ.log(trace[i].toString());
			}
		}
		try{
			//write the Colony information 
			String colFile=outputdir+"/"+imageTitle+"ColonyData.csv";
			File colout = new File(colFile);
			count=0;
			while (colout.exists()){
				count++;
				colFile=outputdir+"/"+imageTitle+"colonyData_"+count+".csv";
				colout= new File(colFile);
			}
			PrintWriter out3 = new PrintWriter(colFile);
			out3.println("Name,Area,Solidity,Perimeter,Circularity,Roundness,Polar Circularity,Min Feret, Max Feret,Minimum Red, Maximum Red, Mean Red, Stdev Red, Skew Red, Kurt Red, 2DcolorSkew Red,Minimum Green, Maximum Green, Mean Green, Stdev Green, Skew Green, Kurt Green, 2DcolorSkew Green,Minimum Blue, Maximum Blue, Mean Blue, Stdev Blue, Skew Blue, Kurt blue, 2DcolorSkew blue, R/G, R/B, G/B,Slice");
			for (int i=0;i<colData.length;i++){
				if (colData[i][1]>0){
					out3.print("Colony "+colData[i][32]+"- Roi"+(i+1));
					for (int j=0;j<colData[0].length;j++){
						out3.print(","+colData[i][j]);	
					}
					out3.println("");
				}
			}
			out3.close();


		}catch (Exception e){
			IJ.log(e.toString());
			StackTraceElement[] trace = e.getStackTrace();
			for (int i=0;i<trace.length;i++){
				IJ.log(trace[i].toString());
			}
		}
	}

	/**
	 * @param ori The Imageplus containing an imageStack of colonies to take the data from
	 * @param outputdir The directory where data is saved
	 */
	public void getColonyHistograms(ImagePlus ori, String outputdir){
		RoiManager rm = RoiManager.getInstance();
		Roi[] crois = rm.getRoisAsArray();
		double [][] histDatar = new double[crois.length][];
		double [][] histDatag = new double[crois.length][];
		double [][] histDatab = new double[crois.length][];
		int [] dim = ori.getDimensions();
		ImagePlus [] histImages = null;
		if (ori.getBitDepth()==16){
			histImages = new ImagePlus[3];
			Duplicator dc = new Duplicator();
			histImages[0]=dc.run(ori,1,1,1,dim[3],1,1);
			histImages[1]=dc.run(ori,2,2,1,dim[3],1,1);
			histImages[2]=dc.run(ori,3,3,1,dim[3],1,1);
			
		} else if (ori.getBitDepth()==24){
			histImages = new ChannelSplitter().split(ori);;
		}	
		for (int k=0;k<histImages.length;k++){
			for (int i=0;i<crois.length;i++){
				ori.setRoi(crois[i]);
				ImageStatistics is= ori.getStatistics();
				if (k==0){
					histDatar[i]=is.histogram();
				} else if (k==1){
					histDatag[i]=is.histogram();
				} else {
					histDatab[i]=is.histogram();
				}
			}
		}
		try {
			//write the pearson correlations
			String rFile=outputdir+"/"+imageTitle+"HistDataRed.csv";
			File rout = new File(rFile);
			String gFile=outputdir+"/"+imageTitle+"HistDataGreen.csv";
			File gout = new File(gFile);
			String bFile=outputdir+"/"+imageTitle+"HistDataBlue.csv";
			File bout = new File(bFile);
			int count=0;
			while (rout.exists()){
				count++;
				rFile=outputdir+"/"+imageTitle+"HistDataRed"+count+".csv";
				rout= new File(rFile);
			}
			count=0;
			while (gout.exists()){
				count++;
				gFile=outputdir+"/"+imageTitle+"HistDataGreen"+count+".csv";
				gout= new File(gFile);
			}
			while (bout.exists()){
				count++;
				bFile=outputdir+"/"+imageTitle+"HistDataBlue"+count+".csv";
				bout= new File(bFile);
			}
			PrintWriter out4 = new PrintWriter(rFile);
			PrintWriter out5 = new PrintWriter(gFile);
			PrintWriter out6 = new PrintWriter(bFile);
			out4.println("\"sep=,\"");
			for (int i=0;i<histDatar.length;i++){
				out4.print("Colony "+(i+1));
				for (int j=0;j<histDatar[0].length;j++){
					out4.print(","+histDatar[i][j]);	
				}
				out4.println("");
			}
			out4.close();

			out5.println("\"sep=,\"");
			for (int i=0;i<histDatag.length;i++){
				out5.print("Colony "+(i+1));
				for (int j=0;j<histDatag[0].length;j++){
					out5.print(","+histDatag[i][j]);	
				}
				out5.println("");
			}
			out5.close();

			out6.println("\"sep=,\"");
			for (int i=0;i<histDatab.length;i++){
				out6.print("Colony "+(i+1));
				for (int j=0;j<histDatab[0].length;j++){
					out6.print(","+histDatab[i][j]);	
				}
				out6.println("");
			}
			out6.close();
			
		} catch (Exception e){
			IJ.log(e.toString());
			StackTraceElement[] trace = e.getStackTrace();
			for (int i=0;i<trace.length;i++){
				IJ.log(trace[i].toString());
			}
		}
	}

	/**
	 * 	
	 * @param l the double list of pearson correlations
	 * @param Cols the amount of columns in the data
	 * @return the string that can be saved that contains all pearson data
	 */
	public String getPearsonTable(java.util.List<Double> l, int Cols){
		String pTable="";
		for (int i=1;i<Cols;i++){
			pTable+=",Colony "+(i+1);
		}
		pTable+="\n";
		int count =0;
		for (int i=0;i<Cols-1;i++){
			pTable+="Colony "+ (i+1);
			for (int k=0;k<i;k++){
				pTable+=",";
			}
			for (int j=i;j<Cols-1;j++){
				pTable+=","+l.get(count);
				count++;
			}
			pTable+="\n";
		}
		return pTable;
	}
	
	/**
	 * 
	 * @param d1 a 2 dimensional float array to compare to d2
	 * @param d2 the 2nd 2-dimensional float array to compare d1 with
	 * @return the pearson correlation coefficient
	 */
	public double getPearson2d (float [][] d1, float [][] d2){	
		int count = (d1.length+1)*d1[0].length+1;
		double sumxy = 0;
		double sumx =0;
		double sumy=0;
		double sumxsq =0;
		double sumysq =0;
		for (int i=0; i<d1.length;i++){
			for (int j=0;j<d1[0].length;j++){
				sumx+=(double)(d1[i][j]);
				sumy+=(double)(d2[i][j]);
				sumxsq+=(double)(d1[i][j]*d1[i][j]);
				sumysq+=(double)(d2[i][j]*d2[i][j]);
				sumxy+=(double)(d1[i][j]*d2[i][j]);
			}
		}
		double p = (count*sumxy-sumx*sumy)/(Math.sqrt(count*sumxsq-sumx*sumx)*Math.sqrt(count*sumysq-sumy*sumy));
		return p;
	}
	
	/**
	 * @param d1 a 1-dimensional float array to correlate with d2
	 * @param d2 the 2nd 1-dimensional float array to compare to d1
	 * @return the pearson correlation coefficient 
	 */
	public double getPearson1D (double [] d1, double [] d2){
		if (d1.length!=d2.length){
			//IJ.log("Arrays are not the same length, the result is incomplete");
		}
		int l = Math.min(d1.length,d2.length);
		int count = (l+1);
		double sumxy = 0;
		double sumx =0;
		double sumy=0;
		double sumxsq =0;
		double sumysq =0;
		
		for (int i=0; i<l;i++){
			sumx+=d1[i];
			sumy+=d2[i];
			sumxsq+=d1[i]*d1[i];
			sumysq+=d2[i]*d2[i];
			sumxy+=d1[i]*d2[i];
		}
		double p = (count*sumxy-sumx*sumy)/(Math.sqrt(count*sumxsq-sumx*sumx)*Math.sqrt(count*sumysq-sumy*sumy));
		return p;
	}

	/**
	 * 
	 * @param im the image of which the maximum difference profile is calculated
	 * @return the double array of the image line that has most variation
	 */
	public double[] getMaxProfile(ImagePlus im){
		double [] Diffs = new double[im.getWidth()];
		double dif = 0;
		IJ.run(im, "Line Width...", "line=5");
		
		for (int i=2; i<im.getHeight()-2;i++){
			double av;
			double sum=0;
			im.setRoi(new Line(0,i,im.getWidth(),i));
			ProfilePlot pf = new ProfilePlot(im);
			double[] data = pf.getProfile();
			for (int j=0;j<data.length;j++){
				sum+=data[j];
			}
			av=sum/data.length;
			double sumdiff=0;
			for (int j=0;j<data.length;j++){
				sumdiff+=Math.abs(data[j]-av);
			}
			sumdiff/=av;
			if (sumdiff>dif){
				Diffs = data;
//				IJ.log(""+i+", "+sumdiff);
				dif=sumdiff;
			}
		}
		return Diffs;
	}
	
	/**
	 * 
	 * @param ori The imagestack on which to calculate the information
	 * @param outputdir The saving directory
	 * @return double array of data that is measured on the colonies
	 */
	public double [][] getColonyInformation(ImagePlus ori, String outputdir){
		ImagePlus OR = ori.duplicate();
		//OR.show();
		double scale=1;
		if (OR.getWidth() >300){
			scale =OR.getWidth()/300;
			IJ.run(OR, "Size...", "width=300 height="+OR.getHeight()/OR.getWidth()*300+" depth="+ori.getNSlices() +" constrain average interpolation=Bilinear");
		}
		ImagePlus [] orisplit=null;
		ImagePlus [] orisplit2=null;
		int [] dim = OR.getDimensions();
		if (OR.getBitDepth()==16){
			IJ.run(OR, "8-bit", "");
			orisplit = new ImagePlus[3];
			Duplicator dc = new Duplicator();
			orisplit[0]=dc.run(OR,1,1,1,dim[3],1,1);
			orisplit[1]=dc.run(OR,2,2,1,dim[3],1,1);
			orisplit[2]=dc.run(OR,3,3,1,dim[3],1,1);
			orisplit2[0]=dc.run(OR,1,1,1,dim[3],1,1);
			orisplit2[1]=dc.run(OR,2,2,1,dim[3],1,1);
			orisplit2[2]=dc.run(OR,3,3,1,dim[3],1,1);

		} else if (OR.getBitDepth()==24){
			orisplit = new ChannelSplitter().split(OR);
			orisplit2 = new ChannelSplitter().split(OR);
			//orisplit[0].show();
			//orisplit[1].show();
			//orisplit[2].show();
			//orisplit2[0].show();
			//orisplit2[1].show();
			//orisplit2[2].show();
		}
		
		int minsize =orisplit[0].getWidth()*orisplit[0].getHeight()/100;
		IJ.run(ori, "Options...", "iterations=3 count=1 black pad do=Nothing");
		
		
		for (int i=0; i<orisplit.length;i++){
			IJ.run(orisplit[i], "Auto Local Threshold", "method=Phansalkar radius="+(OR.getWidth()*2/3)+" parameter_1=0 parameter_2=0 white stack");
			IJ.run(orisplit2[i], "Auto Local Threshold", "method=Contrast radius="+(OR.getWidth()*2/3)+" parameter_1=0 parameter_2=0 white stack");
			IJ.run(orisplit[i], "Analyze Particles...", "size="+minsize+"-Infinity show=Masks  in_situ stack");
			IJ.run(orisplit2[i], "Analyze Particles...", "size="+minsize+"-Infinity show=Masks in_situ stack");//excludeIJ.run(orisplit[i], "Fill Holes", "stack");
			IJ.run(orisplit[i], "Close-", "stack");
			IJ.run(orisplit[i], "Watershed", "stack");
			IJ.run(orisplit2[i], "Fill Holes", "stack");
			IJ.run(orisplit2[i], "Close-", "stack");
			IJ.run(orisplit2[i], "Watershed", "stack");
			IJ.run(orisplit[i], "Analyze Particles...", "size="+minsize+"-Infinity show=Masks  in_situ stack");
			IJ.run(orisplit2[i], "Analyze Particles...", "size="+minsize+"-Infinity show=Masks exclude in_situ stack");//exclude
			
			//orisplit2[i].show();
		}

		ImageCalculator ic = new ImageCalculator();
		ImagePlus i3 =ic.run("AND stack create",orisplit[0],orisplit[1]);
		//i3.show();
		ic = new ImageCalculator();
		ImagePlus i4 =ic.run("AND stack create",i3,orisplit[2]);
		//i4.show();
		i4.setTitle("combined");
		//IJ.run(i4, "Close-", "stack");
		//IJ.run(i4, "Fill Holes", "stack");
		IJ.run(i4, "Analyze Particles...", "size="+minsize+"-Infinity show=Masks in_situ stack");//circularity=0.8-1.00 //exclude

		ImageCalculator ic2 = new ImageCalculator();
		ImagePlus i5 =ic2.run("OR stack create",orisplit2[0],orisplit2[1]);
		//i3.show();
		ic2 = new ImageCalculator();
		ImagePlus i6 =ic2.run("OR stack create",i5,orisplit2[2]);
		//i6.show();
		i6.setTitle("combined2");
		//IJ.run(i6, "Close-", "stack");
		//IJ.run(i6, "Fill Holes", "stack");
		IJ.run(i6, "Analyze Particles...", "size="+minsize+"-Infinity show=Masks in_situ stack");//circularity=0.8-1.00 //exclude
		
		ImageCalculator ic3 = new ImageCalculator();
		ImagePlus i7 =ic3.run("OR stack create",i6,i4);
		IJ.run(i7, "Open", "stack");
		//i7.show();
		i7.setTitle("combined both");
		IJ.run(ori, "Options...", "iterations=1 count=1 black pad do=Nothing");

		// now we have the mask to analyze on the three differen channels
		// now we need to make an roi array from the masks

		RoiManager rm = RoiManager.getInstance();
		if (rm==null){
			rm = new RoiManager();
		} else {
			rm.close();
			rm = new RoiManager();
		}
		//IJ.log(""+minsize);
		IJ.run(i7, "Analyze Particles...", "size="+minsize+"-Infinity show=Masks in_situ add exclude stack");
		Roi[] crois = rm.getRoisAsArray();
		Roi[] convexrois=new Roi[crois.length];
		double[][] outputdata = new double[crois.length][33];
		//new WaitForUserDialog("check").show();
		ArrayList<Integer> removelist = new ArrayList<Integer>();
		for (int j=0;j<crois.length;j++){ // replacing rois with convex hulls
			i7.setRoi(crois[j]);
			outputdata[j][32]=crois[j].getZPosition();
			crois[j]=new PolygonRoi(crois[j].getConvexHull(),Roi.POLYGON);
			crois[j].setPosition(1,1,(int)outputdata[j][32]);
			i7.setRoi(crois[j]);
			ImageStatistics isj = new ImageStatistics();
			i7.setSlice((int)outputdata[j][32]);
			isj= i7.getStatistics(Measurements.AREA+Measurements.ELLIPSE+Measurements.MEAN);
			if (((4*isj.area)/(Math.PI*isj.major*isj.major))>0.7 && (isj.mean/255.0)>0.90 && ((Math.PI*4*isj.area)/(crois[j].getLength()*crois[j].getLength()))>0.89){
				//IJ.log(""+isj.areaFraction+" "+isj.mean);
				convexrois[j]=crois[j];
				//IJ.log(""+convexrois[j].getStatistics().area);
				rm.setRoi(crois[j],j);
			} else {
				//crois[j]=null;
				removelist.add(j);
			}
		}
		if (removelist.size()>0){
			int[] rl = new int[removelist.size()];
			for (int k=0;k<removelist.size();k++){
				rl[k]=removelist.get(k);
			}
			rm.setSelectedIndexes(rl);
			rm.runCommand("Delete");
		}
		//ew WaitForUserDialog("check2").show();
		
		// use the rois to go over the three different colors and get the data we want
		

		IJ.log(""+rm.getCount());
		ImagePlus [] split2 = null;
		dim = ori.getDimensions();
		if (ori.getBitDepth()==16){
			split2 = new ImagePlus[3];
			Duplicator dc = new Duplicator();
			split2[0]=dc.run(ori,1,1,1,dim[3],1,1);
			split2[1]=dc.run(ori,2,2,1,dim[3],1,1);
			split2[2]=dc.run(ori,3,3,1,dim[3],1,1);
			
		} else if (ori.getBitDepth()==24){
			split2 = new ChannelSplitter().split(ori);
		}	
		
		rm.runCommand("Save", outputdir+"/"+imageTitle+"RoiSet.zip");
		for (int i=0; i<split2.length;i++){
			for (int j=0;j<crois.length;j++){
				split2[i].setSlice((int)outputdata[j][32]);
				split2[i].setRoi(convexrois[j]);
				if (convexrois[j]!=null){
					//IJ.log("scale "+scale);
					IJ.run(split2[i], "Scale... ", "x="+scale+" y="+scale);
					ImageStatistics is = new ImageStatistics();
					
					if (i==0){
						is=split2[i].getStatistics(Measurements.AREA+Measurements.ELLIPSE+Measurements.MEAN+Measurements.MIN_MAX+Measurements.STD_DEV+Measurements.KURTOSIS+Measurements.SKEWNESS+Measurements.FERET);
						outputdata[j][0]=convexrois[j].getStatistics().area;
						//Roi cvh = new PolygonRoi(crois[j].getConvexHull(),Roi.POLYLINE); // get the convex hull
						outputdata[j][1]=outputdata[j][0]/convexrois[j].getStatistics().area; // solidity
						//IJ.log(""+convexrois[j].getStatistics().area);
						outputdata[j][2]=convexrois[j].getLength(); // perimeter
						outputdata[j][3]=(4*Math.PI*outputdata[j][0])/(outputdata[j][2]*outputdata[j][2]); // circularity
						outputdata[j][4]=(4*outputdata[j][0])/(Math.PI*crois[j].getStatistics().major*crois[j].getStatistics().major); // roundness
						outputdata[j][5]=new Polar_Circularity().getPolarCircularity(convexrois[j],split2[i]); // polar circularity
						outputdata[j][6]=convexrois[j].getFeretValues()[2]; // min feret
						outputdata[j][7]=convexrois[j].getFeretValues()[0]; // max feret
	
						outputdata[j][8]=is.min;
						outputdata[j][9]=is.max;
						outputdata[j][10]=is.mean;
						outputdata[j][11]=is.stdDev;
						outputdata[j][12]=is.skewness;
						outputdata[j][13]=is.kurtosis;
						double xdis =is.xCentroid-is.xCenterOfMass;
						double ydis = is.yCentroid-is.yCenterOfMass;
						outputdata[j][14]=(Math.sqrt(xdis*xdis+ydis*ydis))/outputdata[j][7];
						
					} else if (i==1) {
						is=convexrois[j].getStatistics();
						outputdata[j][15]=is.min;
						outputdata[j][16]=is.max;
						outputdata[j][17]=is.mean;
						outputdata[j][18]=is.stdDev;
						outputdata[j][19]=is.skewness;
						outputdata[j][20]=is.kurtosis;
						double xdis =is.xCentroid-is.xCenterOfMass;
						double ydis = is.yCentroid-is.yCenterOfMass;
						outputdata[j][21]=(Math.sqrt(xdis*xdis+ydis*ydis))/outputdata[j][7];
					} else {
						is=convexrois[j].getStatistics();
						outputdata[j][22]=is.min;
						outputdata[j][23]=is.max;
						outputdata[j][24]=is.mean;
						outputdata[j][25]=is.stdDev;
						outputdata[j][26]=is.skewness;
						outputdata[j][27]=is.kurtosis;
						double xdis =is.xCentroid-is.xCenterOfMass;
						double ydis = is.yCentroid-is.yCenterOfMass;
						outputdata[j][28]=(Math.sqrt(xdis*xdis+ydis*ydis))/outputdata[j][7];
						outputdata[j][29]=outputdata[j][10]/outputdata[j][17];
						outputdata[j][30]=outputdata[j][10]/outputdata[j][24];
						outputdata[j][31]=outputdata[j][17]/outputdata[j][24];
					}
				}
				
			}
		}
		//IJ.log(Arrays.toString(outputdata[0]));
		return outputdata;
	}
	
	/**
	 * 
	 * @param im The input image from the radial transform
	 * @return an imagePlus that has been filled with the average horizontal value instead of background
	 */
	public ImagePlus getHorizontalAverageFilled(ImagePlus im){
		
		for (int j=0;j<3;j++){
			im.setPosition(j+1);
			ImageProcessor p = im.getChannelProcessor();
			float[][] f = p.getFloatArray();
			//IJ.log(""+f.length);
			for (int i=0;i<f[0].length;i++){ // get the averages per row excluding zeros
				float total =0;
				float count=0;
				for (int h=0;h<f.length;h++){
					if (f[h][i]>0){
						total=total+f[h][i];
						count=count+1;
					}
				}
				float av=total/count;
				for (int h=0;h<f.length;h++){
					if (f[h][i]==0){
						f[h][i]=av;
						//IJ.log(""+i+", "+h+", "+av);
					}
				}
			}
			p.setFloatArray(f);
		}
		return im;
	}
	/**
	 * 
	 * @param im The input image from the radial transform
	 * @return an imagePlus that has been filled with the average vertical value instead of background
	 */
	public ImagePlus getVerticalAverageFilled(ImagePlus im){
		for (int j=0;j<3;j++){
			im.setPosition(j+1);
			ImageProcessor p = im.getChannelProcessor();
			float[][] f = p.getFloatArray();
			//IJ.log(""+f.length);
			for (int i=0;i<f.length;i++){ // get the averages per row excluding zeros
				float total =0;
				float count=0;
				for (int h=0;h<f[i].length;h++){
					if (f[i][h]>0){
						total=total+f[i][h];
						count=count+1;
					}
				}
				float av=total/count;
				for (int h=0;h<f.length;h++){
					if (f[i][h]==0){
						f[i][h]=av;
						//IJ.log(""+i+", "+h+", "+av);
					}
				}
			}
			p.setFloatArray(f);
		}
		return im;
	}
	/**
	 * 
	 * @param pim the original RGB image that has been radially transformed
	 * @return 3 color profiles of the image
	 */
	public ProfilePlot [] getColorProfiles (ImagePlus pim){	
		//IJ.run(pim,"Make Composite", "display=Composite");
		CompositeConverter.makeComposite(pim);
		ProfilePlot [] pfs = new ProfilePlot[3];
		for (int i=0;i<3;i++){
			pim.setPosition(i+1);
			pim.setRoi(new Line(180,0,180,pim.getHeight()));
			pfs[i]=new ProfilePlot(pim);
		}
		return pfs;
	}

	/**
	 * @param imp3 the image to get the 360 degree radial transform of
	 * @return the radial transformed imagePlus
	 */
	public ImagePlus getRadialTransform(ImagePlus imp3){
		int w = imp3.getWidth();
		int h = imp3.getHeight();
		ImagePlus[] channels = null;
		int [] dim = imp3.getDimensions();
		if (imp3.getBitDepth()==16){
			channels = new ImagePlus[3];
			Duplicator dc = new Duplicator();
			channels[0]=dc.run(imp3,1,1,1,dim[3],1,1);
			channels[1]=dc.run(imp3,2,2,1,dim[3],1,1);
			channels[2]=dc.run(imp3,3,3,1,dim[3],1,1);
			
		} else if (imp3.getBitDepth()==24){
			channels = new ChannelSplitter().split(imp3);
		}	

		float [][] red = getRadialKymogram(channels[0],w,h);	
		float [][] green = getRadialKymogram(channels[1],w,h);	
		float [][] blue = getRadialKymogram(channels[2],w,h);	

		ImageStack out = new ImageStack(red.length,red[0].length);
		//IJ.log(""+out.isRGB());
		//IJ.log(""+red.length+" "+red[0].length);
		ImageProcessor ro = new FloatProcessor(red.length,red[0].length);
		ro.setFloatArray(red);
		ImageProcessor go = new FloatProcessor(green);
		ImageProcessor bo = new FloatProcessor(blue);
		out.addSlice("r",ro);
		out.addSlice("g",go);
		out.addSlice("b",bo);
		
		ImagePlus Iout = new ImagePlus("output",out);
		IJ.run(Iout, "8-bit", "");
		
		//Iout.hide();
		return Iout;
	}
	/**
	 * 
	 * @param imp2 the image to get the radial data of
	 * @param w width of the image
	 * @param h height of the image
	 * @return a 2d float array that contains the radial transformed numbers of the image
	 */
	public float[][] getRadialKymogram(ImagePlus imp2,int w,int h){
		imp2.setRoi(new Line(w/2,h/2,w,h/2));
		float [][] vals = new float [360][0];
		int LineL =w/2;
		double [] tmp = new double[w/2];
		for (int i=0;i<360;i++){
			IJ.run(imp2, "Rotate...", "rotate angle=1");
			ProfilePlot profile = new ProfilePlot(imp2);
			tmp = profile.getProfile();	
			float[] tf = new float[LineL];
			for (int j=0;j<Math.min(LineL,tmp.length);j++){
				tf[j]=(float)tmp[j];
			}
			vals[i]=tf;
			//IJ.log(""+tf.length);
		}
		return vals;
		
	}
	
	/**
	 * 
	 * @param imp the imagestack to align
	 * @return the aligned imagestack
	 */
	public ImagePlus [] getAlignedImages(ImagePlus imp){
		ImagePlus OR = imp.duplicate();
		ImagePlus edges = imp.duplicate();
		IJ.run(edges, "Find Edges", "stack");
		IJ.run(imp, "Gaussian Blur...", "sigma=5 stack");
		IJ.run(imp, "8-bit", "");
		IJ.run(imp, "Find Edges", "stack");
		IJ.run(imp, "Gaussian Blur...", "sigma=5 stack");
		IJ.setAutoThreshold(imp, "Moments dark");
		
		IJ.run(imp, "Convert to Mask", "method=Triangle background=Dark calculate black");
		
		for (int i=0; i<imp.getDimensions()[3];i++){ // equalize the images
			imp.setSlice(i+1);
			edges.setSlice(i+1);
			IJ.run(imp, "Enhance Contrast", "saturated=0.01");
			ImageStatistics is = imp.getStatistics(ImageStatistics.CENTER_OF_MASS);
			OR.setSlice(i+1);
			IJ.run(OR, "Translate...", "x="+(imp.getWidth()/2-is.xCenterOfMass)+" y="+(imp.getHeight()/2-is.yCenterOfMass)+" interpolation=None slice");
			IJ.run(edges, "Translate...", "x="+(imp.getWidth()/2-is.xCenterOfMass)+" y="+(imp.getHeight()/2-is.yCenterOfMass)+" interpolation=None slice");
		}
		ImagePlus [] aligned = new ImagePlus[2];
		aligned[0]=OR;
		aligned[1]=edges;
		//imp.hide();
		//imp.close();
		return aligned;
		
	}
	/**
	 * 
	 * @param imp the images from which the colonies are selected
	 * @param wi the imageSize as filled in by the user
	 * @return the stack of selected images, at each point selection a roi of wi by wi has been added to the imagestack
	 */
	public ImagePlus getColonies(ImagePlus imp, int wi){

		IJ.run(imp, "Select None", "");
		IJ.setTool("multipoint");
		WaitForUserDialog wfd = new WaitForUserDialog("Please mark all colonies you want to measure\n\nMake sure the colonies are not touching since this will influence the measurements\npress OK after you are done");
		wfd.show();
		PointRoi pr = new PointRoi(0,0);
		pr = (PointRoi)imp.getRoi();
		while (pr==null){
			wfd = new WaitForUserDialog("You forgot to mark your colonies!\n \nPlease mark all colonies you want to measure\n\nMake sure the colonies are not touching since this will influence the measurements\npress OK after you are done\n \nPress escape to cancel");
			wfd.show();
			pr = (PointRoi)imp.getRoi();
			if (wfd.escPressed()){
				return null;
			}
		}
		imp.hide();
		Point[] ps = pr.getContainedPoints();
		//IJ.log(""+ps.length);
		ImageStack is = new ImageStack(wi,wi);
		ImagePlus b16 = IJ.createHyperStack("colonies",wi,wi,3,ps.length,1,16);
		ImagePlus Icol=null;
		for (int i=0;i<ps.length;i++){
			//IJ.run(imp, "Specify...", "width="+wi+" height="+wi+" x="+ps[i].x+" y="+ps[i].y+" oval constrain centered");
			Roi r = new OvalRoi(ps[i].x-wi/2,ps[i].y-wi/2,wi,wi);
			imp.setRoi(r);
			Duplicator d = new Duplicator();
			if (imp.getBitDepth()==24){
				ImagePlus d2 = d.run(imp);
				is.addSlice(""+i,d2.getProcessor());
			} else if  (imp.getBitDepth()==16){
				ImagePlus d1 = d.run(imp,1,1,1,1,1,1);
				ImagePlus d2 = d.run(imp,2,2,1,1,1,1);
				ImagePlus d3 = d.run(imp,3,3,1,1,1,1);
				b16.setPosition(1,i,1);
				b16.setProcessor(d1.getProcessor());
				b16.setPosition(2,i,1);
				b16.setProcessor(d2.getProcessor());
				b16.setPosition(3,i,1);
				b16.setProcessor(d3.getProcessor());
			}
			//IJ.log(""+wi+", "+ d2.getProcessor().getWidth());
			
		}
		if (imp.getBitDepth()==24){
			Icol = new ImagePlus("colonies",is);
		} else if  (imp.getBitDepth()==16){
			Icol=b16;
		}

		return Icol;
	}
}
