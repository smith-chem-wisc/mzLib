﻿using System;
using System.IO;
using TorchSharp;
using TorchSharp.Modules;

namespace Proteomics.RetentionTimePrediction
{
    /// <summary>
    /// Chronologer is a deep learning model for highly accurate prediction of peptide C18 retention times (reported in % ACN).
    /// Chronologer was trained on a new large harmonized database of > 2.6 million retention time observations
    /// (2.25M unique peptides) constructed from 11 community datasets
    /// and natively supports prediction of 17 different modification types.
    /// With only a few observations of a new modification type (> 10 peptides),
    /// Chronologer can be easily re-trained to predict up to 10 user supplied modifications.
    ///
    /// Damien Beau Wilburn, Ariana E. Shannon, Vic Spicer, Alicia L. Richards, Darien Yeung, Danielle L. Swaney, Oleg V. Krokhin, Brian C. Searle
    /// bioRxiv 2023.05.30.542978; doi: https://doi.org/10.1101/2023.05.30.542978
    /// 
    /// https://github.com/searlelab/chronologer
    ///
    /// Licensed under Apache License 2.0
    /// 
    /// </summary>
    public class Chronologer : torch.nn.Module<torch.Tensor, torch.Tensor>
    {
        public Chronologer() : this(Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
            "RetentionTimePrediction",
            "Chronologer_20220601193755_TorchSharp.dat")) { }

        /// <summary>
        /// Initializes a new instance of the Chronologer model class with pre-trained weights from the paper
        /// Deep learning from harmonized peptide libraries enables retention time prediction of diverse post
        /// translational modifications paper (https://github.com/searlelab/chronologer).
        /// Eval mode is set to true and training mode is set to false by default.
        ///
        /// Please use .Predict() for using the model, not .forward(). 
        /// </summary>
        /// <param name="weightsPath"></param>
        /// <param name="evalMode"></param>
        public Chronologer(string weightsPath, bool evalMode = true) : base(nameof(Chronologer))
        {
            RegisterComponents();

            LoadWeights(weightsPath);//loads weights from the file

            if (evalMode)
            {
                this.eval();
                this.train(false);
            }
        }

        /// <summary>
        /// Do not use for inferring. Use .Predict() instead.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override torch.Tensor forward(torch.Tensor x)
        {
            var input = seq_embed.forward(x).transpose(1, -1);

            var residual = input.clone();
            input = conv_layer_1.forward(input); //renet_block
            input = norm_layer_1.forward(input);
            input = relu.forward(input);
            input = conv_layer_2.forward(input);
            input = norm_layer_2.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);

            residual = input.clone();
            input = conv_layer_4.forward(input);//renet_block
            input = norm_layer_4.forward(input);
            input = relu.forward(input);
            input = conv_layer_5.forward(input);
            input = norm_layer_5.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);

            residual = input.clone();
            input = conv_layer_7.forward(input);//renet_block
            input = norm_layer_7.forward(input);
            input = term_block.forward(input);
            input = relu.forward(input);
            input = conv_layer_8.forward(input);
            input = norm_layer_8.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);
            
            input = dropout.forward(input);
            input = flatten.forward(input);
            input = output.forward(input);
            
            return input;
        }

        /// <summary>
        /// Loads pre-trained weights from the file Chronologer_20220601193755_TorchSharp.dat.
        /// </summary>
        /// <param name="weightsPath"></param>
        private void LoadWeights(string weightsPath)
        {
            //load weights from the file
            this.load(weightsPath, true);
        }

        /// <summary>
        /// Predicts the retention time of the input peptide sequence. The input must be a torch.Tensor of shape (1, 52).
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
        public torch.Tensor Predict(torch.Tensor input)
        {
            return this.call(input);
        }

        //All Modules (shortcut modules are for loading the weights only)
        private Embedding seq_embed = torch.nn.Embedding(55, 64, 0);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_1 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_2 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_3 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_4 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_5 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 2);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_6 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_7 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_8 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 3);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_9 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_1 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_2 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_3 = torch.nn.BatchNorm1d(64); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_4 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_5 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_6 = torch.nn.BatchNorm1d(64); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_7 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_8 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_9 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> term_block = torch.nn.Identity();
        private torch.nn.Module<torch.Tensor, torch.Tensor> relu = torch.nn.ReLU(true);
        private torch.nn.Module<torch.Tensor, torch.Tensor> dropout = torch.nn.Dropout(0.01);
        private torch.nn.Module<torch.Tensor, torch.Tensor> flatten = torch.nn.Flatten(1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> output = torch.nn.Linear(52 * 64, 1);
    }
}
