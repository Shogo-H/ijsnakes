import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.plugin.frame.Recorder;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.TextField;
import java.util.Vector;

import big.ij.snake2D.Snake2DKeeper;
import big.ij.snake2D.Snake2DScale;

/**
 * Implements ExtendedPlugInFilter to use ESnake as a plug-in for ImageJ.
 * 
 * @version November 14, 2013
 * 
 * @author Shogo HIRAMATSU
 */
public class Shonake_ implements ExtendedPlugInFilter {

	/** Image to process. */
	private ImagePlus imp_ = null;

	/** Image processed */
	private ImagePlus dstPlus_ = null;

	/** Initial dialog. */
	private final GenericDialog dialog_ = new GenericDialog("Shonake");

	/** Input image types allowed. */
	private static final int CAPABILITIES = DOES_8G | DOES_16 | DOES_32
			| CONVERT_TO_FLOAT | DOES_STACKS;

	// ----------------------------------------------------------------------------

	/**
	 * Default maximum number of iterations allowed in the snake optimization
	 * process.
	 */
	private static final int DEFAULT_MAX_ITER = 500;

	/** Default number of snake control points. */
	private static final int DEFAULT_NUM_NODES = 4;

	/**
	 * Default standard deviation of the Gaussian prefiltering in the
	 * preprocessing step.
	 */
	private static final int DEFAULT_STD = 3;

	/** Default tradeoff energy parameter. */
	private static final double DEFAULT_ALPHA = 0.00002;

	// ----------------------------------------------------------------------------

	/**
	 * Label for the maximum number of iterations allowed in the snake
	 * optimization process.
	 */
	private static final String MAX_ITER = "Max_iterations";

	/** Label for the number of snake control points. */
	private static final String NUM_NODES = "Control_points";

	/**
	 * Label for the standard deviation of the Gaussian prefiltering in the
	 * preprocessing step.
	 */
	private static final String GAUSSIAN_BLUR = "Gaussian_blur";

	/** Label for the tradeoff energy parameter. */
	private static final String ALPHA = "Alpha";

	/** Label for the stopping criterion of the optimization of the snake. */
	private static final String IMMORTAL = "Immortal";

	/** Label for the brightness type of the object to segment. */
	private static final String TARGET = "Target_brightness";

	/** Label for the energy type to use. */
	private static final String ENERGY_TYPE = "Energy_type";

	/** Label for the saving into the RoiManager of ImageJ. */
	private static final String SAVE = "Save ROI";

	/** List of the brightness types of the object to segment. */
	private static final String[] TARGETOPTIONS = { "Dark", "Bright" };

	/** List of the energy types. */
	private static final String[] ENERGYTYPEOPTIONS = { "Contour", "Region",
			"Mixture" };

	// ----------------------------------------------------------------------------

	/** Maximum number of iterations allowed in the snake optimization process. */
	private static int life_ = DEFAULT_MAX_ITER;

	/** Number of snake control points. */
	private static int M_ = DEFAULT_NUM_NODES;

	/**
	 * Standard deviation of the Gaussian prefiltering in the preprocessing
	 * step.
	 */
	private static int std_ = DEFAULT_STD;

	/** Brightness type of the object to segment.(DETECTDARK or DETECTBRIGHT) */
	private static int detect_ = ESnake.DETECTDARK;

	/** Type of energy used. */
	// private static int energyType_ = ESnake.CONTOURENERGY;
	private static int energyType_ = ESnake.REGIONENERGY;

	/** Tradeoff energy parameter. */
	private static double alpha_ = DEFAULT_ALPHA;

	/** If true, the optimization process continues up to convergence. */
	private static boolean immortalFlag_ = true;

	/** If true, the result is stored in the RoiManager of ImageJ. */
	private static boolean saveROI_ = true;

	/** If true, it is first frame. (made by Shogo HIRAMATSU) */
	private static boolean firstframe_ = true;
	
	/** If 0, it is last frame. (made by Shogo HIRAMATSU) */
	private static int lastframe_ = -1;

	// /** Initial contour. (made by Shogo HIRAMATSU) */
	// private Roi initialContour_ = null;

	/** Stack for processed images. (made by Shogo HIRAMATSU) */
	private ImageStack dstStack_ = null;

	// ============================================================================
	// PUBLIC METHODS

	/**
	 * Filters use this method to process the image. If the SUPPORTS_STACKS flag
	 * was set, it is called for each slice in a stack. With CONVERT_TO_FLOAT,
	 * the filter is called with the image data converted to a FloatProcessor (3
	 * times per image for RGB images). ImageJ will lock the image before
	 * calling this method and unlock it when the filter is finished. For
	 * PlugInFilters specifying the NO_IMAGE_REQUIRED flag and not the DONE
	 * flag, run(ip) is called once with the argument null.
	 */
	@Override
	public void run(ImageProcessor ip) {
		if (firstframe_ == true) {
			/** get numbers in "dialog_" and store to vector from top to bottom */
			final Vector<?> numbers = dialog_.getNumericFields();

			/** get checks in "dialog_" and store to vector from top to bottom */
			final Vector<?> checkboxes = dialog_.getCheckboxes();

			/** get pop-up-menu's choice in "dialog_" and store */
			String targetAsString = dialog_.getNextChoice();
			if (targetAsString.equals(TARGETOPTIONS[ESnake.DETECTDARK])) {
				detect_ = ESnake.DETECTDARK;
			} else if (targetAsString
					.equals(TARGETOPTIONS[ESnake.DETECTBRIGHT])) {
				detect_ = ESnake.DETECTBRIGHT;
			} else {
				IJ.error("Internal error: unexpected brightness detection mode");
				return;
			}

			M_ = (new Integer(((TextField) numbers.elementAt(0)).getText()))
					.intValue();
			std_ = (new Integer(((TextField) numbers.elementAt(1)).getText()))
					.intValue();

			String energytypeAsString = dialog_.getNextChoice();
			if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.CONTOURENERGY])) {
				energyType_ = ESnake.CONTOURENERGY;
			} else if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.REGIONENERGY])) {
				energyType_ = ESnake.REGIONENERGY;
			} else if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.MIXTUREENERGY])) {
				energyType_ = ESnake.MIXTUREENERGY;
			} else {
				IJ.error("Internal error: unexpected energy type selected");
				return;
			}

			alpha_ = (new Double(((TextField) numbers.elementAt(2)).getText()))
					.doubleValue();
			life_ = (new Integer(((TextField) numbers.elementAt(3)).getText()))
					.intValue();

			immortalFlag_ = ((Checkbox) checkboxes.elementAt(0)).getState();
			saveROI_ = ((Checkbox) checkboxes.elementAt(1)).getState();

			// ----------------------------------------------------------------------------

			/** for saving ROI */
			Recorder.setCommand("Shonake ");

			if (targetAsString.equals(TARGETOPTIONS[ESnake.DETECTDARK])) {
				Recorder.recordOption(TARGET, TARGETOPTIONS[ESnake.DETECTDARK]);
			} else if (targetAsString
					.equals(TARGETOPTIONS[ESnake.DETECTBRIGHT])) {
				Recorder.recordOption(TARGET,
						TARGETOPTIONS[ESnake.DETECTBRIGHT]);
			} else {
				IJ.error("Internal error: unexpected brightness detection mode");
				return;
			}

			Recorder.recordOption(NUM_NODES, "" + M_);
			Recorder.recordOption(GAUSSIAN_BLUR, "" + std_);

			if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.REGIONENERGY])) {
				Recorder.recordOption(ENERGY_TYPE,
						ENERGYTYPEOPTIONS[ESnake.REGIONENERGY]);
			} else if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.CONTOURENERGY])) {
				Recorder.recordOption(ENERGY_TYPE,
						ENERGYTYPEOPTIONS[ESnake.CONTOURENERGY]);
			} else if (energytypeAsString
					.equals(ENERGYTYPEOPTIONS[ESnake.MIXTUREENERGY])) {
				Recorder.recordOption(ENERGY_TYPE,
						ENERGYTYPEOPTIONS[ESnake.MIXTUREENERGY]);
			} else {
				IJ.error("Internal error: unexpected energy type selected");
				return;
			}

			Recorder.recordOption(ALPHA, "" + alpha_);
			Recorder.recordOption(MAX_ITER, "" + life_);
			Recorder.recordOption(IMMORTAL, "" + immortalFlag_);
			Recorder.recordOption(SAVE, "" + saveROI_);

			if (saveROI_)
				Recorder.saveCommand();

		}

		// ----------------------------------------------------------------------------

		// ip.setSliceNumber(ip.getSliceNumber()+1);

		ESnake mysnake = new ESnake((FloatProcessor) ip, std_, life_, M_,
				alpha_, immortalFlag_, detect_, energyType_, imp_.getRoi());
		Snake2DKeeper keeper = new Snake2DKeeper();

		// if (IJ.isMacro()) {

		// keeper.optimize(mysnake, imp_);
		keeper.optimize(mysnake, null);

		// } else {
		// keeper.interactAndOptimize(mysnake, imp_);
		// }

		if (!mysnake.isCanceledByUser()) {
			if (saveROI_) {
				// RoiManager roiManager = RoiManager.getInstance();
				// if (roiManager == null){
				// roiManager = new RoiManager();
				// }
				Snake2DScale[] skin = mysnake.getScales();
				PolygonRoi roi = new PolygonRoi(skin[1], Roi.TRACED_ROI);
				// if (saveROI_)
				// roiManager.addRoi(roi);
				ImageProcessor ip2 = ip.duplicate();
				ip2.fill(roi);
				dstStack_.addSlice("", ip2);
				if (firstframe_ == true) {
					dstPlus_ = new ImagePlus("dstStack_", dstStack_);
				}else{
					dstPlus_.setStack(dstStack_);
				}
				dstPlus_.show();
				imp_.setRoi(roi);
			}
		}

		firstframe_ = false;
	}

	// ----------------------------------------------------------------------------

	/**
	 * This method is called by ImageJ to inform the plug-in filter about the
	 * passes to its run method.
	 */
	@Override
	public void setNPasses(int nPasses) {

	}

	// ----------------------------------------------------------------------------

	/**
	 * This method is called once when the filter is loaded. 'arg', which may be
	 * blank, is the argument specified for this plug-in in IJ_Props.txt or in
	 * the plugins.config file of a jar archive containing the plug-in. 'imp' is
	 * the currently active image. This method should return a flag word that
	 * specifies the filters capabilities.
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		imp_ = imp;
		dstStack_ = new ImageStack(imp_.getWidth(), imp_.getHeight());
		return (CAPABILITIES);
	}

	// ----------------------------------------------------------------------------

	/**
	 * This method is called after setup(arg, imp) unless the DONE flag has been
	 * set.
	 */
	@Override
	public int showDialog(final ImagePlus imp, final String command,
			final PlugInFilterRunner pfr) {

		dialog_.addChoice(TARGET, TARGETOPTIONS, TARGETOPTIONS[detect_]);
		dialog_.addNumericField(NUM_NODES, M_, 0);
		dialog_.addNumericField(GAUSSIAN_BLUR, std_, 0);
		dialog_.addChoice(ENERGY_TYPE, ENERGYTYPEOPTIONS,
				ENERGYTYPEOPTIONS[energyType_]);
		dialog_.addNumericField(ALPHA, alpha_, 6);
		dialog_.addNumericField(MAX_ITER, life_, 0);
		dialog_.addCheckbox(IMMORTAL, immortalFlag_);
		dialog_.addCheckbox(SAVE, saveROI_);

		dialog_.addPanel(new ESCreditsButton());

		if (Macro.getOptions() != null) {
			activateMacro(imp);
			return (CAPABILITIES);
		}
		dialog_.showDialog();
		if (dialog_.wasCanceled()) {
			return (DONE);
		}
		if (dialog_.wasOKed()) {
			return (CAPABILITIES);
		} else {
			return (DONE);
		}
	}

	// ============================================================================
	// PRIVATE METHODS

	/**
	 * Prepares the plug-in for running in Macro mode.
	 */
	private void activateMacro(final ImagePlus imp) {
		@SuppressWarnings("unchecked")
		final Vector<TextField> numbers = dialog_.getNumericFields();
		@SuppressWarnings("unchecked")
		final Vector<Choice> choices = dialog_.getChoices();
		@SuppressWarnings("unchecked")
		final Vector<Checkbox> checkboxes = dialog_.getCheckboxes();

		final Choice targetChoice = choices.elementAt(0);
		final TextField numNodes = numbers.elementAt(0);
		final TextField gaussSmooth = numbers.elementAt(1);
		final Choice energytypeChoice = choices.elementAt(1);
		final TextField alphaText = numbers.elementAt(2);
		final TextField iterations = numbers.elementAt(3);
		final Checkbox lifeState = checkboxes.elementAt(0);
		final Checkbox saveState = checkboxes.elementAt(1);
		final String options = Macro.getOptions();

		targetChoice.select(Macro.getValue(options, TARGET, ""
				+ TARGETOPTIONS[detect_]));
		numNodes.setText(Macro.getValue(options, NUM_NODES, "" + M_));
		gaussSmooth.setText(Macro.getValue(options, GAUSSIAN_BLUR, "" + std_));
		energytypeChoice.select(Macro.getValue(options, ENERGY_TYPE, ""
				+ ENERGYTYPEOPTIONS[energyType_]));
		alphaText.setText(Macro.getValue(options, ALPHA, "" + alpha_));
		iterations.setText(Macro.getValue(options, MAX_ITER, "" + life_));

		String s1 = new String(Macro.getValue(options, IMMORTAL, ""
				+ immortalFlag_));
		if (s1.equals("true")) {
			lifeState.setState(true);
		} else {
			lifeState.setState(false);
		}

		String s2 = new String(Macro.getValue(options, SAVE, "" + saveROI_));
		if (s2.equals("true")) {
			saveState.setState(true);
		} else {
			saveState.setState(false);
		}
	}
}
