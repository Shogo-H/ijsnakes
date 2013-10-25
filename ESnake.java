import ij.IJ;
import ij.gui.Roi;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;

import java.awt.Color;
import java.awt.Polygon;
import java.awt.geom.Point2D;
import java.util.*;

import big.ij.snake2D.Snake2D;
import big.ij.snake2D.Snake2DNode;
import big.ij.snake2D.Snake2DScale;

/**
 * Exponential spline snake.
 * 
 * @version June 30, 2011
 * 
 * @author Ricard Delgado-Gonzalo (ricard.delgado@gmail.com)
 */
class ESnake implements Snake2D {

	/** Snake defining nodes. */
	private Snake2DNode[] coef_ = null;

	/** LUT with the samples of the B-spline basis function at rate R. */
	private double[] splineFunc_ = null;
	/**
	 * LUT with the samples of the derivative of the B-spline basis function at
	 * rate R.
	 */
	private double[] splinePrimeFunc_ = null;
	/**
	 * $q(m) = \int_{-\infty}^{\infty}\,phi(t)\,phi(t-m)\,\mathrm{d}t$, where
	 * $\phi$ is the basis function.
	 */
	private double[] q_ = null;
	/** Length of the support of the exponential B-spline basis function */
	private static int N = 3;

	/**
	 * LUT with the samples of the x coordinates of the snake contour at rate R.
	 */
	private double[] xPosSkin_ = null;
	/**
	 * LUT with the samples of the y coordinates of the snake contour at rate R.
	 */
	private double[] yPosSkin_ = null;
	/**
	 * LUT with the samples of the x coordinates of the derivative of the snake
	 * contour at rate R.
	 */
	private double[] xPrimePosSkin_ = null;
	/**
	 * LUT with the samples of the y coordinates of the derivative of the snake
	 * contour at rate R.
	 */
	private double[] yPrimePosSkin_ = null;
	/**
	 * LUT with the samples of the x coordinates of the outer ellipse at rate R.
	 */
	private double[] xPosEllipse_ = null;
	/**
	 * LUT with the samples of the y coordinates of the outer ellipse at rate R.
	 */
	private double[] yPosEllipse_ = null;
	/**
	 * LUT with the samples of the x coordinates of the derivative of the outer
	 * ellipse at rate R.
	 */
	private double[] xPrimePosEllipse_ = null;
	/**
	 * LUT with the samples of the y coordinates of the derivative of the outer
	 * ellipse at rate R.
	 */
	private double[] yPrimePosEllipse_ = null;

	/** Horizontal inferior limit of the bounding box of the snake contour. */
	private int xminS_ = 0;
	/** Horizontal superior limit of the bounding box of the snake contour. */
	private int xmaxS_ = 0;
	/** Vertical inferior limit of the bounding box of the snake contour. */
	private int yminS_ = 0;
	/** Vertical superior limit of the bounding box of the snake contour. */
	private int ymaxS_ = 0;
	/** Horizontal inferior limit of the bounding box of the outer ellipse. */
	private int xminE_ = 0;
	/** Horizontal superior limit of the bounding box of the outer ellipse. */
	private int xmaxE_ = 0;
	/** Vertical inferior limit of the bounding box of the outer ellipse. */
	private int yminE_ = 0;
	/** Vertical superior limit of the bounding box of the outer ellipse. */
	private int ymaxE_ = 0;

	/** LUT with M*R samples of one period of a sine. */
	private double[] sinLUT_ = null;
	/** LUT with M*R samples of one period of a cosine. */
	private double[] cosLUT_ = null;

	/** Horizontal first-order Fourier descriptor of the snake contour. */
	private double[] xFourierS_ = null;
	/** Vertical first-order Fourier descriptor of the snake contour. */
	private double[] yFourierS_ = null;
	/**
	 * Expansion factor dimension-wise applied to the Fourier descriptors to
	 * obtain the outer ellipse.
	 */
	private double lambda_ = 0;

	/**
	 * Weights associated to the sine to associate the snake coefficients to the
	 * first-order horizontal Fourier descriptor.
	 */
	private double[] hs_ = null;
	/**
	 * Weights associated to the cosine to associate the snake coefficients to
	 * the first-order horizontal Fourier descriptor.
	 */
	private double[] hc_ = null;

	/** Signed area of the region enclosed by the snake. */
	private double areaSnake_ = 0;
	/**
	 * Signed area of the region enclosed by elliptical approximation of the
	 * snake.
	 */
	private double areaEllipse_ = 0;

	/** Width of the original image data. */
	private int width_ = 0;
	/** Height of the original image data. */
	private int height_ = 0;
	/** Width of the original image data minus two. */
	private int widthMinusTwo_ = 0;
	/** Height of the original image data minus two. */
	private int heightMinusTwo_ = 0;

	/** Initial contour. */
	private Roi initialContour_ = null;

	/** Original image data. */
	private float[] imageData_ = null;
	/** Original image data filtered with a Laplacian filter. */
	private float[] laplacianImageData_ = null;
	/** Preintegrated image data along the vertical direction. */
	private double[] fuy_ = null;
	/**
	 * Preintegrated and filtered (with a Laplacian filter) image data along the
	 * vertical direction.
	 */
	private double[] fuyLap_ = null;

	/** If true indicates that the snake is able to keep being optimized. */
	private boolean alive_ = true;
	/**
	 * If true indicates that the snake will keep iterating till the optimizer
	 * decides so.
	 */
	private boolean immortal_ = true;
	/**
	 * If true, indicates that the user chose to interactively abort the
	 * processing of the snake. Otherwise, if false, indicates that the dealings
	 * with the snake were terminated without user assistance.
	 */
	private boolean canceledByUser_ = false;

	/** Number of spline vector coefficients. */
	private int M_ = 0;
	/** Energy tradeoff factor. */
	private double alpha_ = 0;
	/** Number of iterations left when the immortalFlag is false. */
	private int life_ = 0;
	/** Maximum number of iterations allowed when the immortalFlag is false. */
	private int maxLife_ = 0;
	/** Indicates the energy function of the snake. **/
	private int energyType_ = 0;

	/** PI. */
	private static double PI = Math.PI;
	/** PI/M. */
	private double PIM_ = 0;
	/** 2*PI/M. */
	private double PI2M_ = 0;
	/** PI*(2*cos(PI/M)/M)^2. */
	private double PI4cos2PIMM2_ = 0;

	/** Sampling rate at which the contours are discretized. */
	private static final int DISCRETIZATIONSAMPLINGRATE = 500;
	/** N*DISCRETIZATIONSAMPLINGRATE. */
	private int NR_ = 0;
	/** M*DISCRETIZATIONSAMPLINGRATE. */
	private int MR_ = 0;

	/** Energy contribution of the inner region. */
	private double internalContribution_ = 0;
	/** Energy contribution of the outer region. */
	private double externalContribution_ = 0;
	/**
	 * LUT for the derivatives of the total energy with respect to the defining
	 * parameters of the snake.
	 */
	private Point2D.Double[] energyGradient_ = null;
	/**
	 * LUT for the derivatives of the contour energy with respect to the
	 * defining parameters of the snake.
	 */
	private Point2D.Double[] contourEnergyGradient_ = null;
	/**
	 * LUT for the derivatives of the region energy with respect to the defining
	 * parameters of the snake.
	 */
	private Point2D.Double[] regionEnergyGradient_ = null;
	/**
	 * LUT for the derivatives of the snake area with respect to the defining
	 * parameters of the snake.
	 */
	private Point2D.Double[] dAs_ = null;
	/**
	 * LUT for the derivatives of the internal contribution with respect to the
	 * defining parameters of the snake.
	 */
	private Point2D.Double[] dInternalContribution_ = null;
	/**
	 * LUT for the derivatives of the external contribution with respect to the
	 * defining parameters of the snake.
	 */
	private Point2D.Double[] dExternalContribution_ = null;
	/**
	 * LUT for the derivatives of lambda with respect to the defining parameters
	 * of the snake.
	 */
	private Point2D.Double[] dLambda_ = null;
	/** LUT with $Q(i,j)=\int_{0}^{M}\,\phi_M(t-j)\,\phi'_M(t-i)\,\mathrm{d}t$. */
	private double[][] tableQ_ = null;

	/** Smallest discretization of a Laplacian filter. */
	private static int LAPLACIAN_KERNEL[] = { 0, -1, 0, -1, 4, -1, 0, -1, 0 };

	/**
	 * Constant that indicates that he brightness of the target to segment is
	 * darker than its surroundings.
	 */
	public static final int DETECTDARK = 0;
	/**
	 * Constant that indicates that he brightness of the target to segment is
	 * brighter than its surroundings.
	 */
	public static final int DETECTBRIGHT = 1;

	/**
	 * Constant that indicates that the contour of the snake is parameterized in
	 * clockwise orientation.
	 */
	public static final int CLOCKWISE = 1;
	/**
	 * Constant that indicates that the contour of the snake is parameterized in
	 * counterclockwise orientation.
	 */
	public static final int COUNTERCLOCKWISE = -1;

	/** Constant that indicates that the energy used is based on contours. */
	public static final int CONTOURENERGY = 0;
	/** Constant that indicates that the energy used is based on regions. */
	public static final int REGIONENERGY = 1;
	/**
	 * Constant that indicates that the energy used is based on a combination of
	 * contours and regions.
	 */
	public static final int MIXTUREENERGY = 2;

	/** Auxiliary Gaussian smoother. */
	private GaussianBlur gaussianBlur_ = new GaussianBlur();

	// ============================================================================
	// PUBLIC METHODS

	/**
	 * Constructor.
	 */
	public ESnake(FloatProcessor in, double sigma, int maxLife, int M,
			double alpha, boolean immortal, int detect, int energytype,
			Roi initialContour) {

		if (M < 3) {
			IJ.error("The minimum number of knots for this basis function is three.");
			return;
		}

		maxLife_ = maxLife;
		M_ = M;
		alpha_ = alpha;
		immortal_ = immortal;
		energyType_ = energytype;
		initialContour_ = initialContour;

		FloatProcessor inputCopy = (FloatProcessor) in.duplicate();
		if (detect == DETECTDARK) {
			inputCopy.invert();
		}

		FloatProcessor lapIm = (FloatProcessor) inputCopy.duplicate();
		gaussianBlur_.blurGaussian(lapIm, sigma, sigma, 0.01);
		lapIm.resetMinAndMax();
		lapIm.convolve3x3(LAPLACIAN_KERNEL);
		laplacianImageData_ = (float[]) lapIm.getPixels();

		NR_ = N * DISCRETIZATIONSAMPLINGRATE;
		MR_ = M * DISCRETIZATIONSAMPLINGRATE;
		PIM_ = Math.PI / (double) M;
		PI2M_ = 2 * PIM_;

		PI4cos2PIMM2_ = 2.0 * Math.cos(Math.PI / ((double) M)) / ((double) M);
		PI4cos2PIMM2_ *= (PI4cos2PIMM2_ * Math.PI);

		width_ = inputCopy.getWidth();
		height_ = inputCopy.getHeight();
		widthMinusTwo_ = width_ - 2;
		heightMinusTwo_ = height_ - 2;

		imageData_ = (float[]) inputCopy.getPixels();

		xPosSkin_ = new double[MR_];
		yPosSkin_ = new double[MR_];
		xPrimePosSkin_ = new double[MR_];
		yPrimePosSkin_ = new double[MR_];

		xPosEllipse_ = new double[MR_];
		yPosEllipse_ = new double[MR_];
		xPrimePosEllipse_ = new double[MR_];
		yPrimePosEllipse_ = new double[MR_];
		xFourierS_ = new double[3];
		yFourierS_ = new double[3];

		energyGradient_ = new Point2D.Double[M];
		contourEnergyGradient_ = new Point2D.Double[M];
		regionEnergyGradient_ = new Point2D.Double[M];
		dAs_ = new Point2D.Double[M];
		dInternalContribution_ = new Point2D.Double[M];
		dExternalContribution_ = new Point2D.Double[M];
		dLambda_ = new Point2D.Double[M];
		tableQ_ = new double[M][M];

		for (int i = 0; i < M; i++) {
			energyGradient_[i] = new Point2D.Double(0.0, 0.0);
			contourEnergyGradient_[i] = new Point2D.Double(0.0, 0.0);
			regionEnergyGradient_[i] = new Point2D.Double(0.0, 0.0);
			dAs_[i] = new Point2D.Double(0.0, 0.0);
			dInternalContribution_[i] = new Point2D.Double(0.0, 0.0);
			dExternalContribution_[i] = new Point2D.Double(0.0, 0.0);
			dLambda_[i] = new Point2D.Double(0.0, 0.0);
		}

		life_ = maxLife_;
		buildLUTs();

		initializeContour();

		computeFourierS();
		updateArea();
		computePosSkin();
	}

	// ----------------------------------------------------------------------------

	/**
	 * The purpose of this method is to compute the energy of the snake.
	 */
	@Override
	public double energy() {

		if (!immortal_) {
			life_--;
			if (life_ == 0)
				alive_ = false;
		}

		double contourEnergy, regionEnergy, Etotal;

		if (xminS_ <= 1 || yminS_ <= 1 || xmaxS_ >= widthMinusTwo_
				|| ymaxS_ >= heightMinusTwo_) {
			Etotal = Double.MAX_VALUE;
		} else {
			if (energyType_ == CONTOURENERGY) {
				contourEnergy = computeContourEnergy();
				Etotal = contourEnergy;
			} else if (energyType_ == REGIONENERGY) {
				if (xminE_ <= 1 || yminE_ <= 1 || xmaxE_ >= widthMinusTwo_
						|| ymaxE_ >= heightMinusTwo_) {
					Etotal = Double.MAX_VALUE;
				} else {
					regionEnergy = computeRegionEnergy();
					Etotal = regionEnergy;
				}
			} else {
				if (xminE_ <= 1 || yminE_ <= 1 || xmaxE_ >= widthMinusTwo_
						|| ymaxE_ >= heightMinusTwo_) {
					Etotal = Double.MAX_VALUE;
				} else {
					contourEnergy = computeContourEnergy();
					regionEnergy = computeRegionEnergy();
					Etotal = alpha_ * contourEnergy + (1 - alpha_)
							* regionEnergy;
				}
			}
		}
		return Etotal;
	}

	// ----------------------------------------------------------------------------

	/**
	 * The purpose of this method is to compute the gradient of the snake energy
	 * with respect to the snake-defining nodes.
	 */
	@Override
	public Point2D.Double[] getEnergyGradient() {

		if (xminS_ <= 1 || yminS_ <= 1 || xmaxS_ >= widthMinusTwo_
				|| ymaxS_ >= heightMinusTwo_) {
			for (int i = 0; i < M_; i++) {
				energyGradient_[i].x = 0.0;
				energyGradient_[i].y = 0.0;
			}
		} else {
			if (energyType_ == CONTOURENERGY) {
				computeContourEnergyGradient();
				return contourEnergyGradient_;
			} else if (energyType_ == REGIONENERGY) {
				if (!(xminE_ <= 1 || yminE_ <= 1 || xmaxE_ >= widthMinusTwo_ || ymaxE_ >= heightMinusTwo_)) {
					computeRegionEnergyGradient();
					return regionEnergyGradient_;
				} else {
					for (int i = 0; i < M_; i++) {
						energyGradient_[i].x = 0.0;
						energyGradient_[i].y = 0.0;
					}
				}
			} else {
				if (!(xminE_ <= 1 || yminE_ <= 1 || xmaxE_ >= widthMinusTwo_ || ymaxE_ >= heightMinusTwo_)) {
					computeContourEnergyGradient();
					computeRegionEnergyGradient();
					for (int i = 0; i < M_; i++) {
						energyGradient_[i].x = alpha_
								* contourEnergyGradient_[i].x + (1.0 - alpha_)
								* regionEnergyGradient_[i].x;
						energyGradient_[i].y = alpha_
								* contourEnergyGradient_[i].y + (1.0 - alpha_)
								* regionEnergyGradient_[i].y;
					}
					return energyGradient_;
				} else {
					for (int i = 0; i < M_; i++) {
						energyGradient_[i].x = 0.0;
						energyGradient_[i].y = 0.0;
					}
				}
			}
		}
		return energyGradient_;
	}

	// ----------------------------------------------------------------------------

	/**
	 * This method provides an accessor to the snake-defining nodes.
	 */
	@Override
	public Snake2DNode[] getNodes() {

		return (coef_);
	}

	// ----------------------------------------------------------------------------

	/**
	 * The purpose of this method is to determine what to draw on screen, given
	 * the current configuration of nodes.
	 */
	@Override
	public Snake2DScale[] getScales() {

		Snake2DScale[] skin;
		if (energyType_ == CONTOURENERGY) {
			skin = new Snake2DScale[2];
		} else {
			skin = new Snake2DScale[3];
		}
		skin[0] = new Snake2DScale(Color.YELLOW, new Color(0, 0, 0, 0), true,
				false);
		skin[1] = new Snake2DScale(Color.RED, new Color(0, 0, 0, 0), true,
				false);

		for (int k = 0; k < M_; k++) {
			skin[0].addPoint((int) Math.round(coef_[k].x),
					(int) Math.round(coef_[k].y));
		}

		int rxt, ryt;

		for (int k = 0; k < MR_; k++) {
			rxt = (int) Math.round(xPosSkin_[k] + 0.5);
			ryt = (int) Math.round(yPosSkin_[k] + 0.5);

			if (rxt < 0) {
				rxt = 0;
			} else if (rxt >= width_) {
				rxt = width_ - 1;
			}

			if (ryt < 0) {
				ryt = 0;
			} else if (ryt >= height_) {
				ryt = height_ - 1;
			}

			skin[1].addPoint(rxt, ryt);
		}
		if (energyType_ != CONTOURENERGY) {
			skin[2] = new Snake2DScale(Color.RED, new Color(0, 0, 0, 0), true,
					false);
			for (int k = 0; k < MR_; k++) {
				rxt = (int) Math.round(xPosEllipse_[k]);
				ryt = (int) Math.round(yPosEllipse_[k]);

				if (rxt < 0) {
					rxt = 0;
				} else if (rxt >= width_) {
					rxt = width_ - 1;
				}

				if (ryt < 0) {
					ryt = 0;
				} else if (ryt >= height_) {
					ryt = height_ - 1;
				}

				skin[2].addPoint(rxt, ryt);
			}
		}
		return (skin);
	}

	// ----------------------------------------------------------------------------

	/**
	 * The purpose of this method is to monitor the status of the snake.
	 */
	@Override
	public boolean isAlive() {

		return alive_;
	}

	// ----------------------------------------------------------------------------

	/**
	 * Sets the status of the snake to alive, and restores the maximum number
	 * iterations to the original one.
	 */
	public void reviveSnake() {

		alive_ = true;
		life_ = maxLife_;
	}

	// ----------------------------------------------------------------------------

	/**
	 * This method provides a mutator to the snake-defining nodes.
	 */
	@Override
	public void setNodes(Snake2DNode[] node) {

		for (int i = 0; i < M_; i++) {
			coef_[i].x = node[i].x;
			coef_[i].y = node[i].y;
		}
		computeFourierS();
		updateArea();
		computePosSkin();
	}

	// ----------------------------------------------------------------------------

	/**
	 * Retrieves the area under the curve determined by the snake.
	 */
	public double getArea() {

		return (Math.abs(areaSnake_));
	}

	// ----------------------------------------------------------------------------

	/**
	 * If true, indicates that the user chose to interactively abort the
	 * processing of the snake. Otherwise, if false, indicates that the dealings
	 * with the snake were terminated without user assistance.
	 */
	public boolean isCanceledByUser() {

		return (canceledByUser_);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Retrieves the orientation of the curve determined by the snake.
	 */
	public int getOrientation() {

		return ((int) Math.signum(areaSnake_));
	}

	// ----------------------------------------------------------------------------

	/**
	 * This method is called when the methods Snake2DKeeper.interact(),
	 * Snake2DKeeper.interactAndOptimize(), and Snake2DKeeper.optimize() are
	 * about to terminate. It provides a report on the current status of this
	 * snake.
	 */
	@Override
	public void updateStatus(boolean canceledByUser, boolean snakeDied,
			boolean optimalSnakeFound, Double energy) {

		canceledByUser_ = canceledByUser;
	}

	// ============================================================================
	// PRIVATE METHODS

	/**
	 * Initializes the snake control points. If the input ImagePlus contains an
	 * area ROI, the method computes the snake control points to fit to the
	 * shape.
	 */
	private void initializeContour() {

		coef_ = new Snake2DNode[M_];

		if (initialContour_ != null) {
			IJ.log(initialContour_.getTypeAsString() + " Roi detected.");
			int type = initialContour_.getType();
			if (type == Roi.RECTANGLE || type == Roi.OVAL
					|| type == Roi.POLYGON || type == Roi.FREEROI
					|| type == Roi.TRACED_ROI) {
				IJ.log("Parsing...");
				Polygon p = initialContour_.getPolygon();
				if (p != null) {
					Point2D.Double[] resampledContour = arcLengthResampling(p,
							M_);
					coef_ = getSplineKnots(resampledContour);
				}
				return;
			} else {
				IJ.log("This type of Roi does not enclose any area.");
			}
		}

		IJ.log("Initializing default shape...");
		int radius = (int) (Math.min((double) width_ / 6, (double) height_ / 6));
		int x0 = width_ / 2;
		int y0 = height_ / 2;

		double K = 2 * (1 - Math.cos(PI2M_))
				/ (Math.cos(PIM_) - Math.cos(3 * PIM_));

		for (int i = 0; i < M_; i++) {
			coef_[i] = new Snake2DNode((int) ((double) x0 + radius * K
					* Math.cos(PIM_ * (2 * i + 3))),
					(int) ((double) y0 + radius * K
							* Math.sin(PIM_ * (2 * i + 3))));
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the Fourier descriptors of the outline of the snake.
	 */
	private void computeFourierS() {

		xFourierS_[0] = 0.0;
		yFourierS_[0] = 0.0;
		xFourierS_[1] = 0.0;
		yFourierS_[1] = 0.0;
		xFourierS_[2] = 0.0;
		yFourierS_[2] = 0.0;

		for (int i = 0; i < M_; i++) {
			xFourierS_[0] += coef_[i].x;
			yFourierS_[0] += coef_[i].y;
		}
		xFourierS_[0] /= (double) M_;
		yFourierS_[0] /= (double) M_;

		for (int i = 0; i < M_; i++) {
			xFourierS_[1] += coef_[i].x * hc_[i];
			xFourierS_[2] += coef_[i].x * hs_[i];
			yFourierS_[1] += coef_[i].y * hc_[i];
			yFourierS_[2] += coef_[i].y * hs_[i];
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * $\int_{0}^{M}\,\varphi_{M}(t-k_2)\,\varphi'_{M}(t-k_1)\,\mathrm{d}t$.
	 */
	private double Q(int k1, int k2) {

		double val = 0;
		for (int i = 0; i < MR_; i++) {

			int index1 = i - k2 * DISCRETIZATIONSAMPLINGRATE;

			while (index1 < 0)
				index1 += MR_;

			while (index1 >= MR_)
				index1 -= MR_;

			int index2 = i - k1 * DISCRETIZATIONSAMPLINGRATE;

			while (index2 < 0)
				index2 += MR_;

			while (index2 >= MR_)
				index2 -= MR_;

			val += ESpline3((double) index1
					/ (double) DISCRETIZATIONSAMPLINGRATE)
					* ESpline3_Prime((double) index2
							/ (double) DISCRETIZATIONSAMPLINGRATE);
		}
		val /= (double) DISCRETIZATIONSAMPLINGRATE;
		return val;
	}

	// ----------------------------------------------------------------------------

	/**
	 * Initializes all LUTs of the class.
	 */
	private void buildLUTs() {

		double currentVal;
		splineFunc_ = new double[NR_];
		splinePrimeFunc_ = new double[NR_];
		fuy_ = new double[width_ * height_];
		fuyLap_ = new double[width_ * height_];
		sinLUT_ = new double[MR_];
		cosLUT_ = new double[MR_];
		hs_ = new double[M_];
		hc_ = new double[M_];

		double PI2MR = PI2M_ / (double) DISCRETIZATIONSAMPLINGRATE;
		for (int i = 0; i < MR_; i++) {
			sinLUT_[i] = Math.sin(PI2MR * i + 1.5);
			cosLUT_[i] = Math.cos(PI2MR * i + 1.5);
		}

		double fuyLap_val, fuy_val;
		for (int x = 0; x < width_; x++) {
			fuy_val = 0.0;
			fuyLap_val = 0;
			for (int y = 0; y < height_; y++) {
				int index = x + width_ * y;
				fuy_val += (double) imageData_[index];
				fuy_[index] = fuy_val;
				fuyLap_val += (double) laplacianImageData_[index];
				fuyLap_[index] = fuyLap_val;
			}
		}

		double tmp = 2 * Math.cos(PIM_) / (double) M_;
		for (int i = 0; i < M_; i++) {
			hs_[i] = tmp * Math.sin(PI2M_ * i);
			hc_[i] = tmp * Math.cos(PI2M_ * i);
		}

		for (int i = 0; i < NR_; i++) {
			currentVal = (double) i / (double) DISCRETIZATIONSAMPLINGRATE;
			splineFunc_[i] = ESpline3(currentVal);
			splinePrimeFunc_[i] = ESpline3_Prime(currentVal);
		}

		int qSize = 2 * N - 1;
		q_ = new double[qSize];
		for (int i = 0; i < qSize; i++) {
			q_[i] = qESplineFunc3(i - N + 1);
		}

		for (int i = 0; i < M_; i++) {
			for (int j = 0; j < M_; j++) {
				tableQ_[i][j] = Q(i, j);
			}
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the contour energy.
	 */
	private double computeContourEnergy() {

		double energy = 0.0;
		double fuyLap_val;
		int x1, x2, y1, y2;
		double DeltaX1, DeltaX2, DeltaY1;

		for (int i = 0; i < MR_; i++) {
			x1 = (int) Math.floor(xPosSkin_[i]);
			y1 = (int) Math.floor(yPosSkin_[i]);

			if (x1 < 1) {
				x1 = 1;
			} else if (x1 > widthMinusTwo_) {
				x1 = widthMinusTwo_;
			}

			if (y1 < 1) {
				y1 = 1;
			} else if (y1 > heightMinusTwo_) {
				y1 = heightMinusTwo_;
			}

			x2 = x1 + 1;
			y2 = y1 + 1;

			DeltaX1 = xPosSkin_[i] - x1;
			DeltaY1 = yPosSkin_[i] - y1;
			DeltaX2 = x2 - xPosSkin_[i];

			fuyLap_val = fuyLap_[x1 + width_ * (y1 - 1)] * DeltaX2
					+ fuyLap_[x2 + width_ * (y1 - 1)] * DeltaX1;
			fuyLap_val += 0.5 * ((laplacianImageData_[x1 + width_ * y1]
					* DeltaX2 + laplacianImageData_[x2 + width_ * y1] * DeltaX1) + (DeltaY1 * (((laplacianImageData_[x1
					+ width_ * y1]
					* DeltaX2 + laplacianImageData_[x2 + width_ * y1] * DeltaX1) * (2 - DeltaY1)) + ((laplacianImageData_[x1
					+ width_ * y2]
					* DeltaX2 + laplacianImageData_[x2 + width_ * y2] * DeltaX1) * DeltaY1))));

			energy += fuyLap_val * xPrimePosSkin_[i];
		}
		energy = energy / ((double) DISCRETIZATIONSAMPLINGRATE)
				* getOrientation();
		return (energy);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the region energy.
	 */
	private double computeRegionEnergy() {

		internalContribution_ = getE(imageData_, fuy_, xPosSkin_, yPosSkin_,
				xPrimePosSkin_);
		externalContribution_ = getE(imageData_, fuy_, xPosEllipse_,
				yPosEllipse_, xPrimePosEllipse_);
		return (((externalContribution_ - internalContribution_) - internalContribution_) / areaSnake_);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes a line integral.
	 */
	private double getE(float[] f, double[] fy, double[] x, double y[],
			double[] xp) {

		double fuy_val;
		int x1, x2, y1, y2;
		double DeltaX1, DeltaX2, DeltaY1, DeltaY2;
		double E = 0.0;
		for (int i = 0; i < MR_; i++) {
			x1 = (int) Math.floor(x[i]);
			y1 = (int) Math.floor(y[i]);

			x2 = x1 + 1;
			y2 = y1 + 1;

			DeltaX1 = x[i] - x1;
			DeltaY1 = y[i] - y1;
			DeltaX2 = x2 - x[i];

			DeltaX1 = x[i] - x1;
			DeltaY1 = y[i] - y1;
			DeltaX2 = 1 - DeltaX1;
			DeltaY2 = 1 - DeltaY1;

			fuy_val = DeltaX2
					* (fy[x1 + width_ * y1] - 0.5 * f[x1 + width_ * y1]
							* DeltaY2 * DeltaY2 + 0.5 * f[x1 + width_ * y2]
							* DeltaY1 * DeltaY1)
					+ DeltaX1
					* (fy[x2 + width_ * y1] - 0.5 * f[x2 + width_ * y1]
							* DeltaY2 * DeltaY2 + 0.5 * f[x2 + width_ * y2]
							* DeltaY1 * DeltaY1);
			E -= fuy_val * xp[i];
		}
		return (E / (double) DISCRETIZATIONSAMPLINGRATE);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the gradient of the gradient energy.
	 */
	private void computeContourEnergyGradient() {

		double gradX, gradY, QfuLap;
		int l_p;
		int orientation = getOrientation();
		for (int k = 0; k < M_; k++) {
			gradX = 0.0;
			gradY = 0.0;
			for (int l = k - N + 1; l <= k + N - 1; l++) {
				l_p = l;
				while (l_p < 0)
					l_p += M_;

				while (l_p >= M_)
					l_p -= M_;

				QfuLap = Q_fu(k, l, laplacianImageData_);
				gradX -= coef_[l_p].y * QfuLap;
				gradY += coef_[l_p].x * QfuLap;
			}
			contourEnergyGradient_[k].x = gradX * orientation;
			contourEnergyGradient_[k].y = gradY * orientation;
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the gradient of the region energy.
	 */
	private void computeRegionEnergyGradient() {

		double E = computeRegionEnergy();
		updatedAs();
		updatedEin();
		updatedEout();
		for (int i = 0; i < M_; i++) {
			regionEnergyGradient_[i].x = (dExternalContribution_[i].x - 2.0
					* dInternalContribution_[i].x - E * dAs_[i].x)
					/ areaSnake_;
			regionEnergyGradient_[i].y = (dExternalContribution_[i].y - 2.0
					* dInternalContribution_[i].y - E * dAs_[i].y)
					/ areaSnake_;
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the contour of the snake from the control points.
	 */
	private void computePosSkin() {

		int index;

		// snake bounding box
		xminS_ = width_ - 1;
		xmaxS_ = 0;
		yminS_ = height_ - 1;
		ymaxS_ = 0;

		// ellipse bounding box
		xminE_ = width_ - 1;
		xmaxE_ = 0;
		yminE_ = height_ - 1;
		ymaxE_ = 0;

		double aux, aux2, xPrimeVal, yPrimeVal, xPosVal, yPosVal;
		for (int i = 0; i < MR_; i++) {
			xPosVal = 0.0;
			yPosVal = 0.0;
			xPrimeVal = 0.0;
			yPrimeVal = 0.0;
			for (int k = 0; k < M_; k++) {
				index = i - k * DISCRETIZATIONSAMPLINGRATE;

				while (index < 0)
					index += MR_;

				while (index >= MR_)
					index -= MR_;

				if (index >= NR_) {
					continue;
				} else {
					aux = splineFunc_[index];
					aux2 = splinePrimeFunc_[index];
				}
				xPosVal += coef_[k].x * aux;
				yPosVal += coef_[k].y * aux;
				xPrimeVal += coef_[k].x * aux2;
				yPrimeVal += coef_[k].y * aux2;
			}
			xPosSkin_[i] = xPosVal;
			yPosSkin_[i] = yPosVal;
			xPosEllipse_[i] = xFourierS_[0] + lambda_
					* (xFourierS_[1] * cosLUT_[i] + xFourierS_[2] * sinLUT_[i]);
			yPosEllipse_[i] = yFourierS_[0] + lambda_
					* (yFourierS_[1] * cosLUT_[i] + yFourierS_[2] * sinLUT_[i]);

			// update the bounding box
			if ((int) Math.floor(xPosSkin_[i]) < xminS_)
				xminS_ = (int) Math.floor(xPosSkin_[i]);
			if ((int) Math.ceil(xPosSkin_[i]) > xmaxS_)
				xmaxS_ = (int) Math.ceil(xPosSkin_[i]);
			if ((int) Math.floor(yPosSkin_[i]) < yminS_)
				yminS_ = (int) Math.floor(yPosSkin_[i]);
			if ((int) Math.ceil(yPosSkin_[i]) > ymaxS_)
				ymaxS_ = (int) Math.ceil(yPosSkin_[i]);

			if ((int) Math.floor(xPosEllipse_[i]) < xminE_)
				xminE_ = (int) Math.floor(xPosEllipse_[i]);
			if ((int) Math.ceil(xPosEllipse_[i]) > xmaxE_)
				xmaxE_ = (int) Math.ceil(xPosEllipse_[i]);
			if ((int) Math.floor(yPosEllipse_[i]) < yminE_)
				yminE_ = (int) Math.floor(yPosEllipse_[i]);
			if ((int) Math.ceil(yPosEllipse_[i]) > ymaxE_)
				ymaxE_ = (int) Math.ceil(yPosEllipse_[i]);

			xPrimePosSkin_[i] = xPrimeVal;
			yPrimePosSkin_[i] = yPrimeVal;
			xPrimePosEllipse_[i] = PI2M_
					* lambda_
					* (-xFourierS_[1] * sinLUT_[i] + xFourierS_[2] * cosLUT_[i]);
			yPrimePosEllipse_[i] = PI2M_
					* lambda_
					* (-yFourierS_[1] * sinLUT_[i] + yFourierS_[2] * cosLUT_[i]);
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the contribution of the inner integral of the region energy.
	 */
	private void updatedEin() {

		double gradX, gradY, Qfu;
		int l_p;
		for (int k = 0; k < M_; k++) {
			gradX = 0.0;
			gradY = 0.0;
			for (int l = k - N + 1; l <= k + N - 1; l++) {
				l_p = l;
				while (l_p < 0)
					l_p += M_;

				while (l_p >= M_)
					l_p -= M_;

				Qfu = Q_fu(k, l, imageData_);
				gradX += coef_[l_p].y * Qfu;
				gradY -= coef_[l_p].x * Qfu;
			}
			dInternalContribution_[k].x = gradX;
			dInternalContribution_[k].y = gradY;
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the contribution of the external integral of the region energy.
	 */
	private void updatedEout() {

		double I = externalContribution_ / (lambda_ * PI2M_);
		updatedLambda();
		for (int i = 0; i < M_; i++) {
			dExternalContribution_[i].x = PI2M_
					* (I * dLambda_[i].x + lambda_ * Qx(i, dLambda_[i].x));
			dExternalContribution_[i].y = -PI2M_
					* (I * dLambda_[i].y - lambda_ * Qy(i, dLambda_[i].y));
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the gradient of the snake area.
	 */
	private void updatedAs() {

		for (int i = 0; i < M_; i++) {
			double EdAsdCx = 0;
			double EdAsdCy = 0;
			for (int k = 0; k < M_; k++) {
				EdAsdCx -= coef_[k].y * tableQ_[i][k];
				EdAsdCy -= coef_[k].x * tableQ_[k][i];
			}
			dAs_[i].x = EdAsdCx;
			dAs_[i].y = EdAsdCy;
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the gradient of the area ratio.
	 */
	private void updatedLambda() {

		double lambda2 = lambda_ * lambda_ / 2.0;
		for (int i = 0; i < M_; i++) {
			double gradX = 0;
			double gradY = 0;
			for (int k = 0; k < M_; k++) {
				gradX += coef_[k].y
						* (tableQ_[i][k] + lambda2 * PI4cos2PIMM2_
								* Math.sin(PI2M_ * (k - i)));
				gradY += coef_[k].x
						* (tableQ_[k][i] + lambda2 * PI4cos2PIMM2_
								* Math.sin(PI2M_ * (i - k)));
			}
			gradX /= (lambda_ * areaEllipse_);
			gradY /= (lambda_ * areaEllipse_);
			dLambda_[i].x = gradX;
			dLambda_[i].y = gradY;
		}
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the area of the curve defined by the snake, the enclosing
	 * ellipse and the area ratio.
	 */
	private void updateArea() {

		double area = 0.0;
		int l_p;
		for (int k = 0; k < M_; k++) {
			int kN = k + N;
			for (int l = k - N + 1; l < kN; l++) {
				l_p = l;
				while (l_p < 0)
					l_p += M_;

				while (l_p >= M_)
					l_p -= M_;

				area += coef_[k].y * coef_[l_p].x * q_[kN - l - 1];
			}
		}
		areaSnake_ = area;
		areaEllipse_ = Math.PI
				* (xFourierS_[1] * yFourierS_[2] - yFourierS_[1]
						* xFourierS_[2]);
		lambda_ = Math.sqrt(2.0 * areaSnake_ / areaEllipse_);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Exponential B-spline of order three.
	 */
	private double ESpline3(double t) {

		double ESplineValue = 0.0;
		double eta = 2 * (1 - Math.cos(PI2M_)) / (PI2M_ * PI2M_);
		if ((t >= 0) & (t <= 1)) {
			ESplineValue = (2 * (1 - Math.cos(PIM_ * t) * Math.cos(PIM_ * t)));
		} else if ((t > 1) & (t <= 2)) {
			ESplineValue = (Math.cos(PI2M_ * (t - 2))
					+ Math.cos(PI2M_ * (t - 1)) - 2 * Math.cos(PI2M_));
		} else if ((t > 2) & (t <= 3)) {
			ESplineValue = (1 - Math.cos(PI2M_ * (t - 3)));
		}
		ESplineValue = ESplineValue / (PI2M_ * PI2M_ * eta);
		return (ESplineValue);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Derivative of the exponential B-spline of order three.
	 */
	private double ESpline3_Prime(double t) {

		double ESplinePrimeValue = 0.0;
		double eta = 2 * (1 - Math.cos(PI2M_)) / (PI2M_ * PI2M_);
		if ((t >= 0) & (t <= 1)) {
			ESplinePrimeValue = Math.sin(PI2M_ * t);
		} else if ((t > 1) & (t <= 2)) {
			ESplinePrimeValue = -(Math.sin(PI2M_ * (t - 2)) + Math.sin(PI2M_
					* (t - 1)));
		} else if ((t > 2) & (t <= 3)) {
			ESplinePrimeValue = Math.sin(PI2M_ * (t - 3));
		}
		ESplinePrimeValue = ESplinePrimeValue / (PI2M_ * eta);
		return (ESplinePrimeValue);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Sampled autocorrelation of the exponential B-spline of order three.
	 */
	private double qESplineFunc3(int l) {

		double value = 0.0;
		switch (l) {
		case -2:
			value = (1.0 / 8.0)
					* (PI * Math.cos(PIM_) * Math.sin(PIM_) - (double) M_ + (double) M_
							* Math.cos(PIM_) * Math.cos(PIM_))
					/ ((double) M_ * (1.0 - 2.0 * Math.cos(PIM_)
							* Math.cos(PIM_) + Math.cos(PIM_) * Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_)));
			break;
		case -1:
			value = -(1.0 / 4.0)
					* (PI * Math.cos(PIM_) * Math.sin(PIM_) + (double) M_ - 3.0
							* (double) M_ * Math.cos(PIM_) * Math.cos(PIM_) + 2.0
							* (double) M_
							* Math.cos(PIM_)
							* Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_))
					/ ((double) M_ * (1.0 - 2.0 * Math.cos(PIM_)
							* Math.cos(PIM_) + Math.cos(PIM_) * Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_)));
			break;
		case 1:
			value = (1.0 / 4.0)
					* (PI * Math.cos(PIM_) * Math.sin(PIM_) + (double) M_ - 3.0
							* (double) M_ * Math.cos(PIM_) * Math.cos(PIM_) + 2.0
							* (double) M_
							* Math.cos(PIM_)
							* Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_))
					/ ((double) M_ * (1.0 - 2.0 * Math.cos(PIM_)
							* Math.cos(PIM_) + Math.cos(PIM_) * Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_)));
			break;
		case 2:
			value = -(1.0 / 8.0)
					* (PI * Math.cos(PIM_) * Math.sin(PIM_) - (double) M_ + (double) M_
							* Math.cos(PIM_) * Math.cos(PIM_))
					/ ((double) M_ * (1.0 - 2.0 * Math.cos(PIM_)
							* Math.cos(PIM_) + Math.cos(PIM_) * Math.cos(PIM_)
							* Math.cos(PIM_) * Math.cos(PIM_)));
			break;
		default:
			value = 0.0;
			break;
		}
		return (value);
	}

	// ----------------------------------------------------------------------------

	/**
	 * $\frac{M}{2\,\pi\,\lambda}\,\left(\frac{\mathrm{d}E_{\mathrm{out}}}{\
	 * mathrm
	 * {d}c_x}-\frac{E_{\mathrm{out}}}{\lambda}\,\frac{\mathrm{d}\lambda}{\
	 * mathrm{d}c_x}\right)$.
	 */
	private double Qx(int i, double dLambdadC) {

		double q = 0.0;
		double tmp1, tmp2, DeltaX1, DeltaY1, DeltaX2, DeltaY2;
		int x1, x2, y1, y2;

		for (int j = 0; j < MR_; j++) {

			x1 = (int) Math.floor(xPosEllipse_[j]);
			y1 = (int) Math.floor(yPosEllipse_[j]);

			x2 = x1 + 1;
			y2 = y1 + 1;

			DeltaX1 = xPosEllipse_[j] - x1;
			DeltaY1 = yPosEllipse_[j] - y1;

			DeltaX2 = x2 - xPosEllipse_[j];
			DeltaY2 = y2 - yPosEllipse_[j];

			tmp1 = -yFourierS_[1] * sinLUT_[j] + yFourierS_[2] * cosLUT_[j];
			tmp2 = (1 / (double) M_
					+ (lambda_ * hc_[i] + xFourierS_[1] * dLambdadC)
					* cosLUT_[j] + (lambda_ * hs_[i] + xFourierS_[2]
					* dLambdadC)
					* sinLUT_[j]);

			q += (imageData_[x1 + width_ * y1] * DeltaX2 * DeltaY2
					+ imageData_[x2 + width_ * y1] * DeltaX1 * DeltaY2
					+ imageData_[x1 + width_ * y2] * DeltaX2 * DeltaY1 + imageData_[x2
					+ width_ * y2]
					* DeltaX1 * DeltaY1)
					* tmp1 * tmp2;
		}
		q = q / (double) DISCRETIZATIONSAMPLINGRATE;
		return (q);
	}

	// ----------------------------------------------------------------------------

	/**
	 * $\frac{M}{2\,\pi\,\lambda}\,\left(\frac{\mathrm{d}E_{\mathrm{out}}}{\
	 * mathrm
	 * {d}c_y}-\frac{E_{\mathrm{out}}}{\lambda}\,\frac{\mathrm{d}\lambda}{\
	 * mathrm{d}c_y}\right)$.
	 */
	private double Qy(int i, double dLambdadC) {

		double q = 0.0;
		double tmp1, tmp2, DeltaX1, DeltaY1, DeltaX2, DeltaY2;
		int x1, x2, y1, y2;

		for (int j = 0; j < MR_; j++) {

			x1 = (int) Math.floor(xPosEllipse_[j]);
			y1 = (int) Math.floor(yPosEllipse_[j]);

			x2 = x1 + 1;
			y2 = y1 + 1;

			DeltaX1 = xPosEllipse_[j] - x1;
			DeltaY1 = yPosEllipse_[j] - y1;

			DeltaX2 = x2 - xPosEllipse_[j];
			DeltaY2 = y2 - yPosEllipse_[j];

			tmp1 = -xFourierS_[1] * sinLUT_[j] + xFourierS_[2] * cosLUT_[j];
			tmp2 = (1 / (double) M_
					+ (lambda_ * hc_[i] + yFourierS_[1] * dLambdadC)
					* cosLUT_[j] + (lambda_ * hs_[i] + yFourierS_[2]
					* dLambdadC)
					* sinLUT_[j]);

			q -= (imageData_[x1 + width_ * y1] * DeltaX2 * DeltaY2
					+ imageData_[x2 + width_ * y1] * DeltaX1 * DeltaY2
					+ imageData_[x1 + width_ * y2] * DeltaX2 * DeltaY1 + imageData_[x2
					+ width_ * y2]
					* DeltaX1 * DeltaY1)
					* tmp1 * tmp2;
		}
		q = q / (double) DISCRETIZATIONSAMPLINGRATE;
		return (q);
	}

	// ----------------------------------------------------------------------------

	/**
	 * $\int_{-\infty}^{\infty}\,f(r(t))\,\varphi(t-k)\,\varphi
	 * '(t-l)\,\mathrm{d}t$.
	 */
	private double Q_fu(int k, int l, float[] array) {

		if (Math.abs(l - k) >= N)
			return (0.0);

		double q_fu = 0.0;
		double tmp1, tmp2, DeltaX1, DeltaY1, DeltaX2, DeltaY2;
		int index, index2;
		int x1, x2, y1, y2;

		for (int i = 0; i < NR_; i++) {
			index = i + (k - l) * DISCRETIZATIONSAMPLINGRATE;
			index2 = i + DISCRETIZATIONSAMPLINGRATE * k;

			while (index2 < 0)
				index2 += MR_;

			while (index2 >= MR_)
				index2 -= MR_;

			x1 = (int) Math.floor(xPosSkin_[index2]);
			y1 = (int) Math.floor(yPosSkin_[index2]);

			x2 = x1 + 1;
			y2 = y1 + 1;

			DeltaX1 = xPosSkin_[index2] - x1;
			DeltaY1 = yPosSkin_[index2] - y1;

			DeltaX2 = x2 - xPosSkin_[index2];
			DeltaY2 = y2 - yPosSkin_[index2];

			tmp1 = splineFunc_[i];

			if (index < 0 || index >= NR_) {
				continue;
			} else {
				tmp2 = splinePrimeFunc_[index];
			}
			q_fu += (array[x1 + width_ * y1] * DeltaX2 * DeltaY2
					+ array[x2 + width_ * y1] * DeltaX1 * DeltaY2
					+ array[x1 + width_ * y2] * DeltaX2 * DeltaY1 + array[x2
					+ width_ * y2]
					* DeltaX1 * DeltaY1)
					* tmp1 * tmp2;
		}
		q_fu = q_fu / (double) DISCRETIZATIONSAMPLINGRATE;
		return (q_fu);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Computes the location of the spline coefficients given an array of points
	 * the spline must interpolate using exponential B-splines of order four.
	 */
	private Snake2DNode[] getSplineKnots(Point2D.Double[] contour) {

		double[] knotsX = new double[M_];
		double[] knotsY = new double[M_];
		double b = ESpline3(1.5);

		for (int i = 0; i < M_; i++) {
			knotsX[i] = contour[i].x;
			knotsY[i] = contour[i].y;
		}

		double[] pole = { (-b + Math.sqrt(2 * b - 1)) / (1 - b) };
		knotsX = prescaledPeriodic(knotsX, pole);
		knotsY = prescaledPeriodic(knotsY, pole);

		Snake2DNode[] newCoeff = new Snake2DNode[M_];
		for (int i = 0; i < M_; i++) {
			newCoeff[i] = new Snake2DNode(knotsX[i], knotsY[i]);
		}
		return newCoeff;
	}

	// ----------------------------------------------------------------------------

	/**
	 * Parameterizes a curve to arc-length parameterization with a given number
	 * of points.
	 */
	private Point2D.Double[] arcLengthResampling(Polygon p, int nPoints) {

		p.addPoint(p.xpoints[0], p.ypoints[0]);

		double[] arcLength = new double[p.npoints];
		arcLength[0] = 0;
		for (int i = 1; i < p.npoints; i++) {
			arcLength[i] = arcLength[i - 1]
					+ Math.sqrt((p.xpoints[i] - p.xpoints[i - 1])
							* (p.xpoints[i] - p.xpoints[i - 1])
							+ (p.ypoints[i] - p.ypoints[i - 1])
							* (p.ypoints[i] - p.ypoints[i - 1]));
		}

		Point2D.Double[] resampledCurve = new Point2D.Double[nPoints];
		double delta = arcLength[p.npoints - 1] / nPoints;
		int index = 0;
		for (int i = 0; i < nPoints; i++) {
			double t = delta * i;
			boolean found = false;
			for (; index < (p.npoints - 1) && !found; index++) {
				if (arcLength[index] <= t && arcLength[index + 1] >= t) {
					found = true;
				}
			}
			index--;
			resampledCurve[i] = new Point2D.Double(
					((arcLength[index + 1] - t) * p.xpoints[index] + (t - arcLength[index])
							* p.xpoints[index + 1])
							/ (arcLength[index + 1] - arcLength[index]),
					((arcLength[index + 1] - t) * p.ypoints[index] + (t - arcLength[index])
							* p.ypoints[index + 1])
							/ (arcLength[index + 1] - arcLength[index]));
		}
		return resampledCurve;
	}

	// ----------------------------------------------------------------------------

	/**
	 * Filters an array with a all-pole recursive filter with periodic boundary
	 * conditions.
	 */
	private double[] prescaledPeriodic(double[] s, double[] pole) {

		final int N = s.length;
		for (int p = 0, P = pole.length; (p < P); p++) {
			final double z = pole[p];
			double z1 = z;
			for (int k = N - 1; (0 < k); k--) {
				s[0] += z1 * s[k];
				z1 *= z;
			}
			s[0] /= 1.0 - z1;
			for (int k = 1; (k < N); k++) {
				s[k] += z * s[k - 1];
			}
			z1 = z;
			final int K = N - 1;
			for (int k = 0; (k < K); k++) {
				s[K] += z1 * s[k];
				z1 *= z;
			}
			s[K] *= 1.0 / (1.0 - z1);
			z1 = 1.0 - z;
			z1 *= z1;
			s[K] *= z1;
			for (int k = N - 2; (0 <= k); k--) {
				s[k] = z * s[k + 1] + z1 * s[k];
			}
		}
		return (s);
	}
	

	// ----------------------------------------------------------------------------

	/**
	 * Detect self-intersection (created by Shogo HIRAMATSU)
	 * If there is self-intersection, this method return true.
	 */
	private boolean selfIntersectionCheck(){
		TreeMap<Double, Integer> Ls = new TreeMap<Double, Integer>();//scanning lines
		for(int i=0; i<MR_; i++){
			Ls.put(xPosSkin_[i], i);
		}
		ArrayList<Integer> Ss = new ArrayList<Integer>();//segments on a scanning line
		Map.Entry<Double, Integer> minl;
		while(Ls.isEmpty()!=false){
			minl = Ls.pollFirstEntry();
			if(minl.getKey() < xPosSkin_[minl.getValue()-1] && minl.getKey() < xPosSkin_[minl.getValue()+1]){
				//誕生
				if(Ss.isEmpty()!=true){
					//何もない所に誕生したとき
					if(yPosSkin_[minl.getValue()-1] > yPosSkin_[minl.getValue()+1]){
						Ss.add(0, minl.getValue());
						Ss.add(1, minl.getValue()-1);
					}else{
						Ss.add(0, minl.getValue()-1);
						Ss.add(1, minl.getValue());
					}
				}else if(yPosSkin_[minl.getValue()] < getYofSonL(Ss.get(0), minl.getKey())){
					//一番上に誕生したとき
					if(yPosSkin_[minl.getValue()-1] > yPosSkin_[minl.getValue()+1]){
						Ss.add(0, minl.getValue());
						Ss.add(1, minl.getValue()-1);
					}else{
						Ss.add(0, minl.getValue()-1);
						Ss.add(1, minl.getValue());
					}
					if(intersectionCheck(Ss.get(1), Ss.get(2))) return true;
				}else if(getYofSonL(Ss.get(Ss.size()-1), minl.getKey()) < yPosSkin_[minl.getValue()]){
					//一番下に誕生したとき
					if(yPosSkin_[minl.getValue()-1] > yPosSkin_[minl.getValue()+1]){
						Ss.add(Ss.size(), minl.getValue());
						Ss.add(Ss.size()+1, minl.getValue()-1);
					}else{
						Ss.add(Ss.size(), minl.getValue()-1);
						Ss.add(Ss.size()+1, minl.getValue());
					}
					if(intersectionCheck(Ss.get(Ss.size()-3), Ss.get(Ss.size()-2))) return true;
				}else{
					//中途半端なところに誕生したとき
					for(int i=0; i<Ss.size()-1; i++){
						if(getYofSonL(Ss.get(i), minl.getKey()) < yPosSkin_[minl.getValue()] && yPosSkin_[minl.getValue()] < getYofSonL(Ss.get(i+1), minl.getKey())){
							if(yPosSkin_[minl.getValue()-1] > yPosSkin_[minl.getValue()+1]){
								Ss.add(i+1, minl.getValue());
								Ss.add(i+2, minl.getValue()-1);
							}else{
								Ss.add(i+1, minl.getValue()-1);
								Ss.add(i+2, minl.getValue());
							}
							if(intersectionCheck(Ss.get(i), Ss.get(i+1))) return true;
							if(intersectionCheck(Ss.get(i+2), Ss.get(i+3))) return true;
						}
						
					}
				}	
			}else if(minl.getKey() > xPosSkin_[minl.getValue()-1] && minl.getKey() > xPosSkin_[minl.getValue()+1]){
				//消滅
				Ss.remove(Ss.indexOf(minl.getValue()-1));
				Ss.remove(Ss.indexOf(minl.getValue()));
				
			}else{
				//誕生も消滅もしない点
			}
				Ls.put(yPosSkin_[minl.getValue()], minl.getValue());
				Ls.put(yPosSkin_[minl.getValue()-1], minl.getValue()-1);

				
					s1 = ABOVE(s, L);
					s2 = BELOW(s, L):
					if	(s1 intersects s) then add (s1, s) to A;
					if	(s2 intersects s) then add (s2, s) to A;
			}
				else if (p is a right endpoint) then
				begin
					s = segment of which p is endpoint;
					s1 = ABOVE(s, L);
					s2 = BELOW(s, L):
					if	(s1 intersects s2 to the right of p) then add (s1, s2) to A;
					DELETE(s, L);
				end
				else  /* p is an intersection */
				begin
					s1, s2 = segments that intersect at p, with s1 above s2 to the left of p.
					s3 = ABOVE(s1 , L);
					s4 = BELOW(s2 , L):
					if	(s3 intersects s2) then add (s3, s2) to A;
					if	(s1 intersects s4) then add (s1, s4) to A;
					Interchange s1 and s2 in L;
				end
				/* The detected intersections must now be processed */
				while (A != phi) do
				begin
					(s1, s2) = get and remove first element of A;
					x = x-coordinate of intersection point of s1 and s2;
					if	(MEMBER(x, E) = FALSE) then
					begin
						report (s1, s2);
						INSERT(x, E);
					end
				end
			end
		return false;
	}
	
	// ----------------------------------------------------------------------------

	/**
	 * Detect intersection (created by Shogo HIRAMATSU)
	 * If 2 sections intersect, this method return true.
	 */
	private boolean intersectionCheck(int s1, int s2){
		double ax, ay, bx, by, cx, cy, dx, dy;
		ax = xPosSkin_[s1];
		ay = yPosSkin_[s1];
		bx = xPosSkin_[s1+1];
		by = yPosSkin_[s1+1];
		cx = xPosSkin_[s2];
		cy = yPosSkin_[s2];
		dx = xPosSkin_[s2+1];
		dy = yPosSkin_[s2+1];
		double caba, bada, bcdc, dcac;
		caba = (cx-ax)*(by-ay)-(cy-ay)*(bx-ax);
		bada = (bx-ax)*(dy-ay)-(by-ay)*(dx-ax);
		bcdc = (bx-cx)*(dy-cy)-(by-cy)*(dx-cx);
		dcac = (bx-cx)*(ay-cy)-(by-cy)*(ax-cx);
		if(caba*bada>0 && bcdc*dcac>0) return true;
		else return false;
	}
	
	// ----------------------------------------------------------------------------

	/**
	 * 
	 * 
	 */
	private double getYofSonL(int s, double x){
		double x0, y0, x1, y1, r;
		x0 = xPosSkin_[s];
		y0 = yPosSkin_[s];
		x1 = xPosSkin_[s+1];
		y1 = yPosSkin_[s+1];
		if(x1!=x0) r = (x-x0)/(x1-x0);
		else r = 0.5;
		return y0*(1.-r)+y1*r;
	}
	
}
