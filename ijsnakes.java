// firstly please activate the window of source picture

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.Polygon;
import ij.plugin.*;
import ij.plugin.frame.*;

import java.lang.Math;

public class My_snake2 implements PlugIn {
	private ImagePlus impl_src = null;
	private int[] inar_lap = null;
	
	private static int ITERATION = 100;
	private static int GAUSSIAN_KARNEL[] = {1,2,1,2,4,2,1,2,1};
	private static int LAPLACIAN_KARNEL[] = {0,-1,0,-1,4,-1,0,-1,0};
	
	public void run(String arg){
		IJ.log("**** Starting ****");
		
		IJ.log("get current imageplus");
		impl_src = WindowManager.getCurrentImage();

		IJ.log("get imageProcessor of the imageplus");
		ImageProcessor impr_src = impl_src.getProcessor();
		ImageProcessor impr_lap = impr_src.duplicate();
		
		IJ.log("convolve");
		impr_lap.convolve3x3(GAUSSIAN_KARNEL);
		impr_lap.convolve3x3(LAPLACIAN_KARNEL);

		/*
		IJ.log("make imageplus for the imageProcessor");
		ImageStack imst_lap = new ImageStack(impr_lap.getWidth(), impr_lap.getHeight());
		imst_lap.addSlice("",impr_lap);
		ImagePlus impl_lap = new ImagePlus("", imst_lap);
		impl_lap.show();
		*/
		
		IJ.log("cast imageProcessor to float array");
		inar_lap = (int[])impr_lap.getPixels();	//int is needed when the picture has RGB

		IJ.log("initialize snake");
		Polygon p;	//java.awt.Polygon
		int xPoly[] = {150, 250, 325, 375, 450};
		int yPoly[] = {150, 100, 125, 225, 350};
		p = new Polygon(xPoly, yPoly, xPoly.length);
		

		IJ.log("Try&Error(4-neighbor)");
		double E_min = computeContourEnergy(p, impr_src);
		double E_test;
		int px_src, py_src, px_dst, py_dst;
		for(int k=0; k<ITERATION; k++){
			for(int n=0; n<p.npoints; n++){
				px_src = p.xpoints[n];
				py_src = p.ypoints[n];
				px_dst = p.xpoints[n];
				py_dst = p.ypoints[n];

				//up
				p.xpoints[n]=px_src;
				p.ypoints[n]=py_src-1;
				E_test = computeContourEnergy(p, impr_src);
				if(E_min > E_test){
					px_dst=px_src;
					py_dst=py_src-1;
					E_min = E_test;
				}

				//right
				p.xpoints[n]=px_src+1;
				p.ypoints[n]=py_src;
				E_test = computeContourEnergy(p, impr_src);
				if(E_min > E_test){
					px_dst=px_src+1;
					py_dst=py_src;
					E_min = E_test;
				}

				//down
				p.xpoints[n]=px_src;
				p.ypoints[n]=py_src+1;
				E_test = computeContourEnergy(p, impr_src);
				if(E_min > E_test){
					px_dst=px_src;
					py_dst=py_src+1;
					E_min = E_test;
				}

				//left
				p.xpoints[n]=px_src-1;
				p.ypoints[n]=py_src;
				E_test = computeContourEnergy(p, impr_src);
				if(E_min > E_test){
					px_dst=px_src-1;
					py_dst=py_src;
					E_min = E_test;
				}

				p.xpoints[n]=px_dst;
				p.ypoints[n]=py_dst;
				IJ.log("E_min = " + E_min);
			}

		}
		PolygonRoi roi_p;
		roi_p = new PolygonRoi(p, Roi.POLYGON);
		impl_src.setRoi(roi_p);
		
		IJ.log("**** Finishing ****");
	}

	private double computeContourEnergy(Polygon p, ImageProcessor impr_src){
		//IJ.log("culcurate the rdot energy");
		//sum of |P(t+1) - P(t)|
		double E_rdot = 0.;	
		for(int i=0; i<p.npoints-1; i++){
			E_rdot += Math.sqrt(Math.pow(p.xpoints[i+1] - p.xpoints[i], 2) + Math.pow(p.ypoints[i+1] - p.ypoints[i], 2));
		}
		E_rdot += Math.sqrt(Math.pow(p.xpoints[0] - p.xpoints[p.npoints-1], 2) + Math.pow(p.ypoints[0] - p.ypoints[p.npoints-1], 2));
		//IJ.log("E_rdot = " + E_rdot);

		//IJ.log("culcurate the r2dot energy");
		//sum of |{P(t+2)-P(t+1)} - {P(t+1)-P(t)}| = |P(t+2)-2*P(t+1)+P(t)|
		double E_r2dot = 0.;
		for(int i=0; i<p.npoints-2; i++){
			E_r2dot += Math.sqrt(Math.pow(p.xpoints[i+2] - 2*p.xpoints[i+1] + p.xpoints[i], 2) + Math.pow(p.ypoints[i+2] - 2*p.ypoints[i+1] + p.ypoints[i], 2));
		}
		E_r2dot += Math.sqrt(Math.pow(p.xpoints[0] - 2*p.xpoints[p.npoints-1] + p.xpoints[p.npoints-2], 2) + Math.pow(p.ypoints[0] - 2*p.ypoints[p.npoints-1] + p.ypoints[p.npoints-2], 2));
		E_r2dot += Math.sqrt(Math.pow(p.xpoints[1] - 2*p.xpoints[0] + p.xpoints[p.npoints-1], 2) + Math.pow(p.ypoints[1] - 2*p.ypoints[0] + p.ypoints[p.npoints-1], 2));
		//IJ.log("E_r2dot = " + E_r2dot);

		//IJ.log("culcurate the image energy");
		// sum of lapracian value
		int E_lapR = 0;
		int E_lapG = 0;
		int E_lapB = 0;
		for (int i=0; i<p.npoints; i++){
			E_lapR += 256 - inar_lap[p.xpoints[i] + impr_src.getWidth() * p.ypoints[i]] & 0xff;
			E_lapG += 256 - (inar_lap[p.xpoints[i] + impr_src.getWidth() * p.ypoints[i]]>>8) & 0xff;
			E_lapB += 256 - (inar_lap[p.xpoints[i] + impr_src.getWidth() * p.ypoints[i]]>>16) & 0xff;
		}
		int E_lap = E_lapR + E_lapG + E_lapB;
		//IJ.log("E_lap = " + E_lap);

		//IJ.log("culcurate the total energy including energy of constraint");
		double E_total = E_rdot + E_r2dot + (double)E_lap;
		for (int i=0; i<p.npoints; i++){
			if (p.xpoints[i]<=2 || p.ypoints[i]<=2 || p.xpoints[i]>=impr_src.getWidth()-2 || p.ypoints[i]>=impr_src.getHeight()-2){
				E_total = Double.MAX_VALUE;
			}
		}
		//IJ.log("E_total = " + E_total);
		return E_total;
	}

}
