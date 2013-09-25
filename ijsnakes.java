// firstly please activate the window of source picture

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.Polygon;
import ij.plugin.*;
import ij.plugin.frame.*;

import java.lang.Math;

public class My_snake implements PlugIn {
	private ImagePlus impl_src = null;
	private int[] inar_lap = null;
	
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

		IJ.log("make imageplus for the imageProcessor");
		ImageStack imst_lap = new ImageStack(impr_lap.getWidth(), impr_lap.getHeight());
		imst_lap.addSlice("",impr_lap);
		ImagePlus impl_lap = new ImagePlus("", imst_lap);
		impl_lap.show();

		IJ.log("cast imageProcessor to float array");
		inar_lap = (int[])impr_lap.getPixels();	//int is needed when the picture has RGB
		
		IJ.log("set polygon roi");
		Polygon p;	//java.awt.Polygon
		int xPoly[] = {150, 250, 325, 375, 450};
		int yPoly[] = {150, 100, 125, 225, 350};
		p = new Polygon(xPoly, yPoly, xPoly.length);
		Roi po_roi = new PolygonRoi(p, Roi.POLYGON);
		impl_lap.setRoi(po_roi);

		IJ.log("culcurate the vdot energy");
		double E_vdot = 0.;	
		for(int i=0; i<p.npoints-1; i++){
			E_vdot += Math.sqrt(Math.pow(p.xpoints[i+1] - p.xpoints[i], 2) + Math.pow(p.ypoints[i+1] - p.ypoints[i], 2));
		}
		E_vdot += Math.sqrt(Math.pow(p.xpoints[0] - p.xpoints[p.npoints-1], 2) + Math.pow(p.ypoints[0] - p.ypoints[p.npoints-1], 2));
		IJ.log("E_vdot = " + E_vdot);

		IJ.log("culcurate the v2dot energy");
		double E_v2dot = 0.;
		
		
		IJ.log("**** Finishing ****");
	}



}
