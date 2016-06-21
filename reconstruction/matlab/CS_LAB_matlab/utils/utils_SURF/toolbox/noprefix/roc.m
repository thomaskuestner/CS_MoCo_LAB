function varargout = roc(varargin)
%VL_ROC   ROC curve.
%   [TPR,TNR] = VL_ROC(LABELS, SCORES) computes the Receiver Operating
%   Characteristic (ROC) curve. LABELS are the ground truth labels,
%   greather than zero for a positive sample and smaller than zero for
%   a negative one. SCORES are the scores of the samples obtained from
%   a classifier, where lager scores should correspond to positive
%   labels.
%
%   Samples are ranked by decreasing scores, starting from rank 1.
%   TPR(K) and TNR(K) are the true positive and true negative rates
%   when samples of rank smaller or equal to K-1 are predicted to be
%   positive. So for example TPR(3) is the true positive rate when the
%   two samples with largest score are predicted to be
%   positive. Similarly, TPR(1) is the true positive rate when no
%   samples are predicted to be positive, i.e. the constant 0.
%
%   Set the zero the lables of samples that should be ignored in the
%   evaluation. Set to -INF the scores of samples which are not
%   retrieved. If there are samples with -INF score, then the ROC curve
%   may have maximum TPR and TNR smaller than 1.
%
%   [TPR,TNR,INFO] = VL_ROC(...) returns an additional structure INFO
%   with the following fields:
%
%   info.auc:: Area under the ROC curve (AUC).
%     The ROC curve has a `staircase shape' because for each sample
%     only TP or TN changes, but not both at the same time. Therefore
%     there is no approximation involved in the computation of the
%     area.
%
%   info.eer:: Equal error rate (EER).
%     The equal error rate is the value of FPR (or FNR) when the ROC
%     curves intersects the line connecting (0,0) to (1,1).
%
%   VL_ROC(...) with no output arguments plots the ROC curve in the
%   current axis.
%
%   VL_ROC() acccepts the following options:
%
%   Plot:: []
%     Setting this option turns on plotting unconditionally. The
%     following plot variants are supported:
%
%     tntp:: Plot TPR against TNR (standard ROC plot).
%     tptn:: Plot TNR against TPR (recall on the horizontal axis).
%     fptp:: Plot TPR against FPR.
%     fpfn:: Plot FNR against FPR (similar to DET curve).
%
%   NumPositives:: []
%   NumNegatives:: []
%     If set to a number, pretend that LABELS contains this may
%     positive/negative labels. NUMPOSITIVES/NUMNEGATIVES cannot be
%     smaller than the actual number of positive/negative entrires in
%     LABELS. The additional positive/negative labels are appended to
%     the end of the sequence, as if they had -INF scores (not
%     retrieved). This is useful to evaluate large retrieval systems in
%     which one stores ony a handful of top results for efficiency
%     reasons.
%
%   About the ROC curve::
%     Consider a classifier that predicts as positive all samples whose
%     score is not smaller than a threshold S. The ROC curve represents
%     the performance of such classifier as the threshold S is
%     changed. Formally, define
%
%       P = overall num. of positive samples,
%       N = overall num. of negative samples,
%
%     and for each threshold S
%
%       TP(S) = num. of samples that are correctly classified as positive,
%       TN(S) = num. of samples that are correctly classified as negative,
%       FP(S) = num. of samples that are incorrectly classified as positive,
%       FN(S) = num. of samples that are incorrectly classified as negative.
%
%     Consider also the rates:
%
%       TPR = TP(S) / P,      FNR = FN(S) / P,
%       TNR = TN(S) / N,      FPR = FP(S) / N,
%
%     and notice that by definition
%
%       P = TP(S) + FN(S) ,    N = TN(S) + FP(S),
%       1 = TPR(S) + FNR(S),   1 = TNR(S) + FPR(S).
%
%     The ROC curve is the parametric curve (TPR(S), TNR(S)) obtained
%     as the classifier threshold S is varied in the reals. The TPR is
%     also known as recall (see VL_PR()).
%
%     The ROC curve is contained in the square with vertices (0,0) The
%     (average) ROC curve of a random classifier is a line which
%     connects (1,0) and (0,1).
%
%     The ROC curve is independent of the prior probability of the
%     labels (i.e. of P/(P+N) and N/(P+N)).
%
%   REFERENCES:
%   [1] http://en.wikipedia.org/wiki/Receiver_operating_characteristic
%
%   See also: VL_PR(), VL_DET(), VL_HELP().
[varargout{1:nargout}] = vl_roc(varargin{:});
