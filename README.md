Download Link: https://assignmentchef.com/product/solved-stats215-assignment-2
<br>






<strong>Problem 1: </strong><em>Bernoulli GLMs as a latent variable models.</em>

Consider a Bernoulli regression model,

<em>w </em>∼ ~ (<em>µ</em>, <em>Σ</em>)

<em>y<sub>n</sub></em> | <em>x<sub>n</sub></em>,<em>w </em>∼ Bern(<em>f </em>(<em>w</em><sup>T</sup><em>x<sub>n</sub></em>)) for <em>n </em>= 1,…,<em>N</em>,

where <em>w</em>and <em>x<sub>n</sub></em> are vectors in R<em><sup>D</sup></em>, <em>y<sub>n</sub></em> E {0, 1}, and <em>f </em>: R → [0, 1] is the mean function. In class we studied Newton’s method for finding the maximum a posteriori (MAP) estimate <em>w<sup>*</sup></em> = arg max <em>p</em>(<em>w </em>| {<em>x<sub>n</sub></em>, <em>y</em><em>n</em>}<em>N n</em>=1). Now we will consider methods for approximating the full posterior distribution.

<ul>

 <li>Rather than using the logistic function, let the mean function be the normal cumulative distribution function (CDF), or “probit” function,</li>

</ul>

<em>f </em>(<em>u</em>) = Pr(<em>z </em>≤ <em>u</em>) where <em>z </em>∼ iV (0, 1)

=    ~ (<em>z</em>; 0, 1) d<em>z</em>.

~ −∞

<em>u</em>

This is called the probit regression model. Show that the likelihood <em>p</em>(<em>y<sub>n</sub></em> | <em>x<sub>n</sub></em>, <em>w</em>) is a marginal of a joint distribution,

<em>p</em>(<em>y<sub>n</sub></em>,<em>z<sub>n</sub></em> | <em>x<sub>n</sub></em>,<em>w</em>) = 1I[<em>z<sub>n</sub></em> ≥ 0]~[<em>y</em><em>n</em>=1] 1I[<em>z<sub>n</sub></em><em> &lt; </em>0]~<sup>[</sup><em><sup>yn</sup></em><sup>=</sup><sup>0</sup><sup>]</sup>JV (<em>z<sub>n</sub></em> | <em>x</em>T <em>n </em><em>w</em>,1).

&lt;Your answer here.&gt;&gt;

<ul>

 <li>Derive the conditional distributions <em>p</em>(<em>w </em>| {<em>x<sub>n</sub></em>, <em>y</em><em>n</em>,<em>z</em><em>n</em>}<em>N <sub>n</sub></em><sub>=</sub><sub>1</sub>) and <em>p</em>(<em>z<sub>n</sub></em> | <em>x<sub>n</sub></em>, <em>y<sub>n</sub></em>,<em>w</em>).<sup>1</sup></li>

</ul>

&lt;Your answer here.&gt;&gt;

<ul>

 <li><em>Gibbs sampling </em>is a Markov chain Monte Carlo (MCMC) method for approximate posterior inference. It works by repeatedly sampling from the conditional distribution of one variable, holding all others fixed. For the probit regression model, this means iteratively performing these two steps:</li>

</ul>

<ol>

 <li>Sample <em>z<sub>n</sub></em> ∼ <em>p</em>(<em>z<sub>n</sub></em> | <em>x<sub>n</sub></em>, <em>y<sub>n</sub></em>,<em>w</em>) for <em>n </em>= 1,.. .,<em>N </em>holding <em>w </em>fixed;</li>

 <li>Sample <em>w </em>∼ <em>p</em>(<em>w </em>| {<em>x<sub>n</sub></em>, <em>y<sub>n</sub></em>, <em>z</em><em>n</em>}<em>N <sub>n</sub></em><sub>=</sub><sub>1</sub>) holding {<em>z</em><em>n</em>}<em>N <sub>n</sub></em><sub>=</sub><sub>1</sub></li>

</ol>

<sup>1</sup>Observe that <em>z</em><em><sub>n</sub></em> is conditionally independent of {<em>x</em><em><sub>n</sub></em>~, <em>y</em><em><sub>n</sub></em>~, <em>z</em><em>n</em>~}<em>n</em>‘~=<em>n </em>given <em>w</em>.




Note the similarity to EM: rather than computing a posterior distribution over<em> z<sub>n</sub></em>, we draw a sample from it; rather than setting<em> w</em> to maximize the ELBO, we draw a sample from its conditional distribution. It can be shown that this algorithm defines a Markov chain on the space of (<em>w</em>, {<em>z</em><em>n</em>}<em>Nn</em>=1) whose stationary distribution is the posterior<em> p</em>(<em>w</em>, {<em>z</em><em>n</em>}<em>Nn</em>=1 | {<em>x</em><em>n</em>,<em><sub> yn</sub></em><sub>}</sub><em><sub>Nn</sub></em><sub>=</sub><sub>1</sub><sub>)</sub><sub>.</sub> In other words, repeating these steps infinitely many times would yield samples of<em> w</em> and<sub> {</sub><em><sub>zn</sub></em><sub>}</sub><em><sub>Nn</sub></em><sub>=</sub><sub>1</sub> drawn from their posterior distribution.

Implement this Gibbs sampling algorithm and test it on a synthetic dataset with<em> D</em> = 2 dimensional covariates and<em> N</em> = 100 data points. Scatter plot your samples of<em> w</em> and, for comparison, plot the true value of<em> w</em> that generated the data. Do your samples look approximately Gaussian distributed? How does the posterior distribution change when you vary<em> N</em>?

«Your figures and captions here.»

(d)<strong> Bonus.</strong> There are also auxiliary variable methods for logistic regression, where<em> f </em>(<em>u</em>) =<em> e<sup>u</sup></em><em>/</em>(1 +<em><sub> e</sub></em><em>u</em>). Specifically, we have that,

<table>

 <tbody>

  <tr>

   <td width="57"><em><sub>e</sub></em><em>y</em><em>n</em>·<em>w</em><sup>T</sup><em>x</em><em>n</em></td>

   <td rowspan="2" width="455">= f     <em>n</em> exp {( <em>y</em> − 1 <em>x</em>T <em>n </em><em>w</em>−<sup> 1</sup><sub>2</sub><em>z<sub>n</sub></em>(<em>w</em><sup>T</sup><em>x<sub>n</sub></em>)<sup>2</sup>} PG(<em>z<sub>n</sub></em>; 1, 0)d<em>z<sub>n</sub></em>,20</td>

  </tr>

  <tr>

   <td width="57">1 + <em><sub>e</sub></em><em>w</em>T<em> x<sub>n</sub></em></td>

  </tr>

 </tbody>

</table>




where PG(<em>z</em>;<em> b</em>,<em> c</em>) is the density function of the<em> Pólya-gamma</em> (PG) distribution over<em> z</em> ∈ R<sub>+</sub> with parameters<em> b</em> and<em> c</em>. The PG distribution has a number of nice properties: it is closed under exponential tilting so that,

<em><sub>e</sub></em>− 2 1<em>zc</em>2 PG(<em>z</em>;<em> b</em>,0) ∝ PG(<em>z</em>;<em> b</em>,<em> c</em>), and its expectation is available in closed form,

2<em>c </em>tanh ( 2).

Use these properties to derive an EM algorithm for finding<em> w</em> = argmax<em> p</em>({<em>y<sub>n</sub></em>} | {<em>x<sub>n</sub></em>},<em>w</em>). How do the EM updates compare to Newton’s method?




<strong>Problem 2:</strong><em> Spike sorting with mixture models</em>

As discussed in class, “spike sorting” is ultimately a mixture modeling problem. Here we will study the problem in more detail. Let<sub> {</sub><em><sub>y</sub></em><em><sub>n</sub></em><sub>}</sub><em><sub>Nn</sub></em><sub>=</sub><sub>1</sub> represent a collection of spikes. Each<em> y<sub>n</sub></em> ∈ R<em>D</em> is a vector containing features of the<em> n</em>-th spike waveform. For example, the features may be projections of the spike waveform onto the top<em> D</em> principal components. We have the following, general model,

<em>z<sub>n</sub></em> |<em> π </em>∼ <em>π</em>

<em>y<sub>n</sub></em> |<em> z<sub>n</sub></em>,<em> θ </em>∼ <em>p</em>(<em>y<sub>n</sub></em> | <em>θ</em><em>z<sub>n</sub></em>).

The label<em> z<sub>n</sub></em> ∈ {1,…,<em> K</em>} indicates which of the<em> K</em> neurons generated the<em> n</em>-th spike waveform. The probability vector<em> π</em> ∈<em> Δ</em><em>K</em> specifies a prior distribution on spike labels, and the parameters<em> θ</em> = {<em>θ</em><em>k</em>}<em>K k</em>=1 determine the likelihood of the spike waveforms<em> y<sub>n</sub></em> for each of the<em> K</em> neurons. The goal is to infer a posterior distribution<em> p</em>(<em>z<sub>n</sub></em> |<em> y<sub>n</sub></em>,<em> π</em>,<em> θ</em>) over labels for each observed spike, and to learn the parameters<em> π</em><em>~</em> and<em> θ<sup>*</sup></em> that maximize the likelihood of the data.

<ul>

 <li>Start with a Gaussian observation model,</li>

</ul>

<em>y<sub>n</sub></em> |<em> z<sub>n</sub></em>,<em>θ </em>∼ * (<em>y<sub>n</sub></em> | <em>µ</em><em>z<sub>n</sub></em>, <em>Σ</em><em>z<sub>n</sub></em>),

where<em> θ</em><em>k</em> = (<em>µ</em><em>k</em>, <em>Σ</em><em>k</em>) includes the mean and covariance for the<em> k</em>-th neuron.

Derive an EM algorithm to compute<em> π</em><em>t</em>,<em> θ<sup>*</sup></em> = argmax <em>p</em>({<em>y</em><em>n</em>}<em><sup>N</sup></em><em><sub>n</sub></em>=1 |<em> π</em>,<em> θ</em>). Start by deriving the “responsibilities”<em> w</em><em>nk</em> =<em> p</em>(<em>z<sub>n</sub></em> =<em> k</em> |<em> y<sub>n</sub></em>,<em> π</em><sup>‘</sup>,<em> θ</em><sup>‘</sup>) for fixed parameters<em> π</em><sup>‘</sup> and<em> θ</em><sup>‘</sup>. Then use the responsibilities to compute the expected log joint probability,

2(<em>π</em>,<em>θ</em>) = ~<em>N </em>E<em>p</em>(<em>z</em><em>n</em>|<em>y</em><em>n</em>,<em>π</em>/,<em>θ</em>/) [log <em>p</em>(<em>y<sub>n</sub></em>,<em>z<sub>n</sub></em> |<em> π</em>,<em> θ</em>)].

<em>n</em>=1

Finally, find closed-form expressions for<em> π</em><em>~</em> and<em> θ<sup>*</sup></em> that optimize 2 (<em>π</em>,<em> θ</em>). «Your answer here.»

<ul>

 <li>The Gaussian model can be sensitive to outliers and lead spikes from one neuron to be split into two One way to side-step this issue is to replace the Gaussian with a heavier-tailed distribution like the multivariate Student’s t, which has probability density,</li>

</ul>

<em>Γ</em><u>[(</u><em><u>α</u></em><u>0 </u><u>+ </u><em><u>D</u></em><u>)</u><em><u>/</u></em><u>2</u><u>] </u><u>            1 </u>

<em>p</em>(<em>y<sub>n</sub></em> | <em>θ</em><em>z<sub>n</sub></em>) =                  1<em>/</em>2 I<em>y</em><em>n </em>− <em>µ</em><em>z<sub>n</sub>        </em><em>z</em><em>n </em><em>y<sub>n</sub></em> − <em>µ<sub>z</sub></em>

<em>Γ</em>(<em>α</em>0<em>/</em>2)<em>α<sup>/</sup></em><sup>2</sup><em>π<sup>D</sup></em><em>/</em><sup>2</sup>         L                                                                   .

<em>α</em>0

We will treat<em> α</em>0 as a fixed hyperparameter.

Like the negative binomial distribution studied in HW1, the multivariate Student’s t can also be represented as an infinite mixture,

<em>p</em>(<em>y<sub>n</sub></em> | <em>θ</em>) = J<em> p</em>(<em>y<sub>n</sub></em>,<em>τ<sub>n</sub></em> | <em>θ</em>) d<em>τ<sub>n</sub></em> = J X(<em>y</em><em>n</em>;<em>µ</em><em>z</em><em>n</em>,<em>τ</em>−<em>n</em>1<em>Σ</em><em>z</em><em>n</em>)Gamma(<em>τ</em><em>n</em>; <em>α</em>02 , <sup>1</sup><sub> 2</sub>) d<em>τ<sub>n</sub></em>. We will derive an EM algorithm to find<em> π</em><em>t</em>,<em> θ<sup>*</sup></em> in this model.




First, show that the posterior takes the form

<em>p</em>(<em>τ<sub>n</sub></em>,<em>z<sub>n</sub></em> |<em> y<sub>n</sub></em>,<em>π</em>,<em>θ</em>) =<em> p</em>(<em>z<sub>n</sub></em> |<em> y<sub>n</sub></em>,<em>π</em>,<em>θ</em>) <em>p</em>(<em>τ<sub>n</sub></em> | <em>z<sub>n</sub></em>, <em>y<sub>n</sub></em>, <em>θ</em>)

<em>K                                                             </em>1[[<em>z</em>=<em>k</em>]

HI <em>w</em><em>nk </em>Gamma(<em>τ</em> | <em>a</em><em>nk</em>, <em>b</em><em>nk</em>)]

<em>k</em>

and solve for the parameters<em> w</em><em>nk</em>,<em> a</em><em>nk</em>,<em> b</em><em>nk</em> in terms of<em> y<sub>n</sub></em>,<em> π</em>, and<em> θ</em>. &lt;Your answer here.&gt;&gt;

<ul>

 <li>Now compute the expected log joint probability,</li>

</ul>

~ (<em>π</em>, <em>θ</em>) =  E<em>p</em>(<em>τ</em><em>n</em>,<em>z</em><em>n</em>|<em>y</em><em>n</em>,<em>π</em>~,<em>θ</em>~)[log <em>p</em>(<em>y<sub>n</sub></em>,<em>z<sub>n</sub></em>,<em> τ<sub>n</sub></em> |<em> π</em>,<em> θ</em>)],

<em>n</em>=1

using the fact that E[<em>X</em>] =<em> a/b</em> for<em> X </em>∼ Gamma(<em>a</em>,<em> b</em>). You may omit terms that are constant with respect to<em> π</em> and<em> θ</em>.

&lt;Your answer here.&gt;&gt;

<ul>

 <li>Finally, solve for<em> π</em><em>~</em> and<em> θ<sup>*</sup></em> that maximize the expected log joint probability. How does your answer compare to the solution you found in part (a)?</li>

</ul>

&lt;Your answer here.&gt;&gt;




<strong>Problem 3:</strong><em> Poisson matrix factorization</em>

Many biological datasets come in the form of matrices of non-negative counts. RNA sequencing data, neural spike trains, and network data (where each entry indicate the number of connections between a pair of nodes) are all good examples. It is common to model these counts as a function of some latent features of the corresponding row and column. Here we consider one such model, which decomposes a count matrix into a superposition of non-negative row and column factors.

Let<em> Y</em> ∈ N<em><sup>M</sup></em><sup>×</sup><em><sup>N</sup></em> denote an observed<em> M</em> ×<em> N</em> matrix of non-negative count data. We model this matrix as a function of non-negative row factors<em> U</em> ∈ R<em><sup>M</sup></em><sup>×</sup><em><sup>K</sup></em><sub> +</sub>and column factors<em> V</em> ∈ R<em>N</em>×<em>K</em>

+ . Let<em> u<sub>m</sub></em> ∈ R<em><sup>K</sup></em><sub> +</sub>and<em> v<sub>n</sub></em> ∈ ~<em>K </em>


denote the<em> m</em>-th and<em> n</em>-th rows of<em> U</em> and<em> V</em>, respectively. We assume that each observed count<em> y</em><em>mn</em> is conditionally independent of the others given its corresponding row and column factors. Moreover, we assume a linear Poisson model,

<em>y</em><em>mn</em> |<em> u<sub>m</sub></em>, <em>v<sub>n</sub></em> ∼ Poisson(<em>u</em><sup>T</sup><em><sub>m</sub></em><em>v<sub>n</sub></em>).

(Since<em> u<sub>m</sub></em> and<em> v<sub>n</sub></em> are non-negative, the mean parameter is valid.) Finally, assume gamma priors,

<em>u</em><em>mk </em>∼ Gamma(<em>α</em>0,<em>β</em>0),

<em>v</em><em>nk </em>∼ Gamma(<em>α</em>0,<em>β</em>0).

Note that even though the gamma distribution is conjugate to the Poisson, here we have an inner product of two gamma vectors producing one Poisson random variable. The posterior distribution is more complicated. The entries of<em> u<sub>m</sub></em> are not independent under the posterior due to the “explaining away” effect. Nevertheless, we will derive a mean-field variational inference algorithm to approximate the posterior distribution.

<ul>

 <li>First we will use an augmentation trick based on the additivity of Poisson random variables; i.e. the fact that</li>

</ul>

<em>y </em>∼ Poisson ( <em>λ</em>)⇐⇒<em>y</em>= <em>y</em><em>k </em>where <em>y</em><em>k </em>∼ Poisson<em>k</em>) independently, <em>k</em>

for any collection of non-negative rates<em> λ</em>1, . . . ,<em>λ</em><em>K</em> ∈ R<sub>+</sub>. Use this fact to write the likelihood <em>p</em>(<em>y</em><em>mn</em> |<em> u<sub>m</sub></em>,<em> v<sub>n</sub></em>) as a marginal of<sub> ajoint</sub> distribution<em> p</em>(<em>y<sub>mn</sub></em>,<sup> ¯</sup><em>y</em><em>mn</em> |<em> u<sub>m</sub></em>,<em> v<sub>n</sub></em>) where<sup> ¯</sup><em>y</em><em>mn</em> = (<em>y</em><em>mn</em>1,… ,<em> y</em><em>mnK</em>) is a length-<em>K</em> vector of non-negative counts. (Hint: this is similar to Problem 1 in that<em> y</em><em>mn</em> is deterministic given ¯<em>y<sub>mn</sub></em>.)

«Your answer here.»

<ul>

 <li>Let<em> Y</em>¯ ∈ N<em>M</em>×<em>N</em>×<em>K</em> denote the augmented data matrix with entries<em><sub> y</sub></em><em><sub>mnk</sub></em> as above. We will use mean field variational inference to approximate the posterior as,</li>

</ul>

<table>

 <tbody>

  <tr>

   <td width="81"><em>p</em>(</td>

   <td width="101"><em>Y</em>¯ , <em>U</em>, <em>V </em>| <em>Y </em>) ≈ <em>q</em>(</td>

   <td width="451"><em>Y</em><sup>¯</sup>)<em> q</em>(<em>U</em>)<em> q</em>(<em>V</em>) = 1 ii<em>M</em> ii<em>q</em>(<em>y<sub>mn</sub></em>) ] 1 <em>m</em>=1ii<em>M</em> ii<em><sup>K</sup></em><em>q</em>(<em>u<sub>m</sub></em><em>k</em>) ] 1 <em>n</em>=1ii ii<em>K</em><em>q</em>(<em>v</em><em>nk</em>)].<em>m</em>=1<em> n</em>=1                                                                 <em>k</em>=1                             <em>k</em>=1</td>

  </tr>

 </tbody>

</table>




We will solve for the optimal posterior approximation via coordinate descent on the KL divergence to the true posterior. Recall that holding all factors except for<em> q</em>(¯<em>y<sub>mn</sub></em>) fixed, the KL is minimized when

<em>q</em>(¯<em>y</em>) ∝ exp { E<em>q</em><em>mn</em>)<em>q</em>(<em>U</em>)<em>q</em>(<em>V</em>) [log <em>p</em>(<em>Y</em>, <em>Y</em><sup>¯</sup>,<em> U</em>,<em> V</em>)]},




where <em>q</em>(I<em>‘<sub>n</sub></em>) = <em>fl</em>(<em>m</em>‘,<em>n</em>‘)(<em>m</em>,<em>n</em>) <em>q</em>(¯<em>y<sub>m</sub></em>‘<em><sub>n</sub></em>‘) denotes all variational factors except for the (<em>m</em>,<em> n</em>)-th. Show that the optimal<em> q</em>(¯<em>y<sub>mn</sub></em>) is a multinomial of the form,

<em>q</em>(¯<em>y<sub>mn</sub></em>) = Mult(¯<em>y<sub>mn</sub></em>;<em> y</em><em>mn</em>,<em> π</em><em>mn</em>),

and solve for<em> π</em><em>mn</em> E <em>Δ</em><em>K</em>. You should write your answer in terms of expectations with respect to the other variational factors.

&lt;Your answer here.&gt;&gt;

<ul>

 <li>Holding all factors but<em> q</em>(<em>u<sub>m</sub></em><em>k</em>) fixed, show that optimal distribution is</li>

</ul>

<em>q</em>(<em>u<sub>m</sub></em><em>k</em>) = Gamma(<em>u</em><em>mk</em>;<em> α</em><em>mk</em>,<em> β</em><em>mk</em>).

Solve for<em> α</em><em>mk</em>,<em> β</em><em>mk</em>; write your answer in terms of expectations with respect to<em> q</em>(¯<em>y<sub>mn</sub></em>) and<em> q</em>(<em>v<sub>n</sub></em><em>k</em>). &lt;Your answer here.&gt;&gt;

<ul>

 <li>Use the symmetry of the model to determine the parameters of the optimal gamma distribution for<em> q</em>(<em>v<sub>n</sub></em><em>k</em>), holding<em> q</em>(¯<em>y<sub>mn</sub></em>) and<em> q</em>(<em>u<sub>m</sub></em><em>k</em>) fixed,</li>

</ul>

<em>q</em>(<em>v<sub>n</sub></em><em>k</em>) = Gamma(<em>v<sub>n</sub></em><em>k</em>;<em> α</em><em>nk</em>,<em> β</em><em>nk</em>).

Solve for<em> α</em><em>nk</em>,<em> β</em><em>nk</em>; write your answer in terms of expectations with respect to<em> q</em>(¯<em>y<sub>mn</sub></em>) and<em> q</em>(<em>u<sub>m</sub></em><em>k</em>). &lt;Your answer here.&gt;&gt;

<ul>

 <li>Now that the form of all variational factors has been determined, compute the required expectations (in closed form) to write the coordinate descent updates in terms of the other variational parameters. Use the fact that E[log<em>X</em>] =<em> ψ</em>(<em>α</em>)log<em> β</em> for <em>X</em> Gamma(<em>α</em>,<em> β</em>), where <em>ψ </em>is the digamma function.</li>

</ul>

&lt;Your answer here.&gt;&gt;

<ul>

 <li>Suppose that<em> Y</em> is a sparse matrix with only<em> S</em> &lt; <em>MN</em> non-zero entries. What is the complexity of this mean-field coordinate descent algorithm?</li>

</ul>

&lt;Your answer here.&gt;&gt;




<strong>Problem 4: </strong><em>Apply Poisson matrix factorization to C. elegans connectomics data</em>

Make a copy of this Colab notebook: <a href="https://colab.research.google.com/drive/1ZMwcB6vzVaXz4WJiNT514b7zB5s3_SBk">https://colab.research.google.com/drive/1ZMwcB6vzVaXz4WJiNT514b7zB5s3_SBk</a>

Use your solutions from Problem 3 to finish the incomplete code cells. Once you’re done, run all the code cells, save the notebook in .ipynb format, print a copy in .pdf format, and submit these files along with the rest of your written assignment.