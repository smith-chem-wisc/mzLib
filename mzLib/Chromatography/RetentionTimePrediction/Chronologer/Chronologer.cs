using System.Reflection;
using TorchSharp;
using TorchSharp.Modules;

namespace Chromatography.RetentionTimePrediction.Chronologer;

/// <summary>
/// Chronologer deep learning model for retention time prediction.
/// Wraps TorchSharp neural network with pre-trained weights from the paper
/// "Deep learning from harmonized peptide libraries enables retention time prediction of diverse post
/// translational modifications" (https://github.com/searlelab/chronologer). 
/// </summary>
internal sealed class Chronologer : torch.nn.Module<torch.Tensor, torch.Tensor>, IDisposable
{
    private bool _disposed;

    /// <summary>
    /// Initializes a new instance of the Chronologer model class with pre-trained weights.
    /// Eval mode is set to true and training mode is set to false by default.
    /// Please use .Predict() for inference, not .forward() directly.
    /// </summary>
    /// <param name="weightsPath">Path to the pre-trained weights file</param>
    /// <param name="evalMode">If true, sets model to evaluation mode (default: true)</param>
    public Chronologer(string? weightsPath = null, bool evalMode = true)
        : base(nameof(Chronologer))
    {
        if (weightsPath == null)
        {
            // Locate embedded weights file
            var info = Assembly.GetExecutingAssembly().GetName();
            using var resourceStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{info.Name}.RetentionTimePrediction.Chronologer.Chronologer_20220601193755_TorchSharp.dat");
            
            if (resourceStream == null)
            {
                weightsPath = Path.Combine(AppContext.BaseDirectory, "Resources", "Chronologer_20220601193755_TorchSharp.dat");

                if (!File.Exists(weightsPath))
                {
                    throw new FileNotFoundException($"Default Chronologer weights file not found at expected location: {weightsPath}. " +
                        $"Please ensure the resource file is correctly embedded or provide a valid path to the weights file.");
                }
            }
            else
            {
                // Extract embedded resource to temporary file
                weightsPath = Path.Combine(Path.GetTempPath(), "Chronologer_20220601193755_TorchSharp.dat");
                using var fileStream = new FileStream(weightsPath, FileMode.Create, FileAccess.Write);
                resourceStream.CopyTo(fileStream);
            }
        }

        RegisterComponents();
        LoadWeights(weightsPath);

        if (evalMode)
        {
            eval(); // Evaluation mode doesn't update weights
            train(false);
        }
    }

    /// <summary>
    /// Do not use for inferring. Use .Predict() instead. 
    /// Why forward() is not used when predicting outside the training method? 
    /// -> https://stackoverflow.com/questions/58508190/in-pytorch-what-is-the-difference-between-forward-and-an-ordinary-method
    /// </summary>
    public override torch.Tensor forward(torch.Tensor x)
    {
        var input = seq_embed.forward(x).transpose(1, -1);

        // First residual block
        var residual = input.clone();
        input = conv_layer_1.forward(input);
        input = norm_layer_1.forward(input);
        input = relu.forward(input);
        input = conv_layer_2.forward(input);
        input = norm_layer_2.forward(input);
        input = relu.forward(input);
        input = term_block.forward(input);
        input = residual + input;
        input = relu.forward(input);

        // Second residual block
        residual = input.clone();
        input = conv_layer_4.forward(input);
        input = norm_layer_4.forward(input);
        input = relu.forward(input);
        input = conv_layer_5.forward(input);
        input = norm_layer_5.forward(input);
        input = relu.forward(input);
        input = term_block.forward(input);
        input = residual + input;
        input = relu.forward(input);

        // Third residual block
        residual = input.clone();
        input = conv_layer_7.forward(input);
        input = norm_layer_7.forward(input);
        input = term_block.forward(input);
        input = relu.forward(input);
        input = conv_layer_8.forward(input);
        input = norm_layer_8.forward(input);
        input = relu.forward(input);
        input = term_block.forward(input);
        input = residual + input;
        input = relu.forward(input);

        // Output layers
        input = dropout.forward(input);
        input = flatten.forward(input);
        input = output.forward(input);

        return input;
    }

    /// <summary>
    /// Loads pre-trained weights from the file.
    /// </summary>
    private void LoadWeights(string weightsPath)
    {
        load(weightsPath, strict: true);
    }

    /// <summary>
    /// Predicts the retention time of the input peptide sequence. 
    /// The input must be a torch.Tensor of shape (1, 52).
    /// </summary>
    internal torch.Tensor Predict(torch.Tensor input)
    {
        if (_disposed)
            throw new ObjectDisposedException(nameof(Chronologer));

        return call(input);
    }

    // All Modules (shortcut modules are for loading the weights only, not used but required for the weights)
    private Embedding seq_embed = torch.nn.Embedding(55, 64, 0);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_1 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_2 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 1);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_3 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1); // shortcut
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_4 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_5 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 2);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_6 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2); // shortcut
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_7 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_8 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 3);
    private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_9 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3); // shortcut
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_1 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_2 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_3 = torch.nn.BatchNorm1d(64); // shortcut
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_4 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_5 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_6 = torch.nn.BatchNorm1d(64); // shortcut
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_7 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_8 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_9 = torch.nn.BatchNorm1d(64);
    private torch.nn.Module<torch.Tensor, torch.Tensor> term_block = torch.nn.Identity();
    private torch.nn.Module<torch.Tensor, torch.Tensor> relu = torch.nn.ReLU(inplace: true);
    private torch.nn.Module<torch.Tensor, torch.Tensor> dropout = torch.nn.Dropout(0.01);
    private torch.nn.Module<torch.Tensor, torch.Tensor> flatten = torch.nn.Flatten(1);
    private torch.nn.Module<torch.Tensor, torch.Tensor> output = torch.nn.Linear(52 * 64, 1);

    protected override void Dispose(bool disposing)
    {
        if (!_disposed)
        {
            _disposed = true;
        }
        base.Dispose(disposing);
    }
}
