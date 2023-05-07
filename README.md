Download Link: https://assignmentchef.com/product/solved-math567-homework-4
<br>
. <em>SOR method </em>The FD2-Poisson handout posted on the course webpage contains Matlab code for numerically solving the Poisson equation

∇<sup>2</sup><em>u </em>= <em>f</em>(<em>x,y</em>)<em>,            </em>(<em>x,y</em>) ∈ Ω = (<em>a,b</em>) × (<em>a,b</em>) <em>u</em>(<em>x,y</em>) = <em>g</em>(<em>x,y</em>) (<em>x,y</em>) ∈ <em>∂</em>Ω <em>,</em>

with a second order accurate FD method. The function fd2poisson from this handout sets up and solves the resulting linear system using full dense Gaussian elimination (you will fix this in the next problem).

Your goal is to implement a similar function called fd2poissonsor in Matlab or whatever language you choose that instead numerically solves the Poisson equation using successive over relaxation (SOR) with the optimal relaxation parameter of <em>ω </em>= 2<em>/</em>(1 + sin(<em>πh</em>)). Your code should produce a solution with a relative residual of ≤ 10<sup>−8</sup>. This means that the residual, which is given as

<em>r</em><em>j,k </em>= −4<em>u</em><em>j,k </em>+ <em>u</em><em>j</em>−1<em>,k </em>+ <em>u</em><em>j</em>+1<em>,k </em>+ <em>u</em><em>j,k</em>−1 + <em>u</em><em>j,k</em>+1 − <em>h</em>2<em>f</em><em>j,k,</em>

must satisfy

<em>.</em>

Additionally, your code should not form any matrices involving the 5-point FD stencil. You may want to refer to Sections 4.1 &amp; 4.2 of the book for examples of how the Gauss-Seidel and SOR methods can be implemented in Matlab.

Illustrate that your code works by solving the Poisson equation from the FD2-Poisson handout with <em>m </em>= 2<sup>7 </sup>−1. Produce a plot of the numerically computed solution and of the error in the solution (i.e. difference between computed solution and true solution as done in the code given in the handout for the direct solver). Include the the plots and your fd2poissonsor code in your homework write-up.

<ol start="2">

 <li><em>Comparing Poisson solvers </em>In this problem you will compare various methods for solving the associated linear system for the second-order accurate finite-difference (FD2) approximation to the Poisson equation:</li>

</ol>

∇<sup>2</sup><em>u</em>(<em>x,y</em>) = 10<em>π</em><sup>2</sup><em>e</em><sup>sin(2<em>π</em>(<em>x</em>+2<em>y</em>))</sup>(−2sin(2<em>π</em>(<em>x </em>+ 2<em>y</em>)) + cos(4<em>π</em>(<em>x </em>+ 2<em>y</em>)) + 1)

(1) <em>u</em>(<em>x,y</em>) = <em>e</em>sin(2<em>π</em>(<em>x</em>+2<em>y</em>)) for (<em>x,y</em>) ∈ <em>∂</em>Ω<em>.</em>

This is the problem that is setup in the Matlab file testFdPoisson on the course webpage.

<ul>

 <li>Download the Matlab function m from the course webpage and modify the code so that it uses Matlab’s <em>sparse matrix </em>capabilities (alternatively, you may implement this function in whatever language you choose, but it must use sparse matrix libraries). This means you should create sparse versions of the D2x and D2y matrices in this function instead of dense versions as the code currently does. Note that to get full credit you should <em>not </em>call the functions sparse or gallery(’poisson’,m,m) in Matlab, instead you should figure out how to construct the D2x and D2y matrices yourself using spdiags and kron (and possibly toeplitz), if using Matlab. Save your modified function as fd2poissonsp. Turn in a listing of your code.</li>

 <li>Download the Matlab functions fd2poissondst and fd2poissonmg from the course webpage (alternatively, you may implement these function in whatever language you choose). These functions solve the FD2 linear system for the Poisson equation using the (fast) discrete sine transform (DST) and classical multigrid (using a damped Jacobi smoother and a “v-cycle”), respectively. Note that you will also need to download the updated dst and idst files from the course webpage for this code to work.</li>

 <li>Modify the testFdPoisson function so that it solves the Poisson equation (1) using:

  <ul>

   <li>the standard dense Gaussian elimination solver (fd2poisson),</li>

   <li>the SOR solver from problem 1 (fd2poissonsor),</li>

   <li>the sparse Gaussian elimination solver (fd2poissonsp),</li>

   <li>the DST based solver (fd2poissondst), and</li>

   <li>the multigrid solver (fd2poissonmg).</li>

  </ul></li>

</ul>

Time how long each method takes to solve the Poisson equation for <em>m </em>= 2<em><sup>k </sup></em>− 1, <em>k </em>= 4<em>,</em>5<em>,…,</em>10. Produce a table that shows the timing results for each method and for each value of <em>m</em>. Note that for the dense solver you will not be able to go much above <em>k </em>= 7 without possibly causing serious issues for your machine. Try to estimate the timings in this case based on the previous values of <em>m</em>. Also, you should run all the timing comparisons a few times to capture the mean behavior of the timing results. Report the specification of the machine you used to do the computations (i.e. make, processor type and speed, and memory).

<ul>

 <li>Which method appears to be the best in terms of actual wall-clock time and in terms of the rate of increase with <em>m</em>?</li>

</ul>

Note that the wall-clock time for the MG method are higher than you would see if the method was implemented in a compiled language and if a better smoother was used (e.g. Red/Black Gauss-Seidel) [1]. Also, the DST method is not as efficient as it could be since it uses an FFT of twice the needed length to compute the solutions. Furthermore, this code uses the DST in both the <em>x </em>and <em>y </em>directions to solve the problem. More efficient methods exist by only transforming in one direction and then solving several de-coupled tridiagonal systems [2].

<ol start="3">

 <li><em>Fast Poisson solver with Neumann boundary conditions </em>Consider Poisson’s equation with zero Neumann boundary conditions:</li>

</ol>

∇<sup>2</sup><em>u </em>= <em>f</em>(<em>x,y</em>)<em>,            </em>(<em>x,y</em>) ∈ Ω = (<em>a,b</em>) × (<em>a,b</em>) <strong>n </strong>· ∇<em>u</em>(<em>x,y</em>) = 0<em>,               </em>(<em>x,y</em>) ∈ <em>∂</em>Ω <em>,</em>

where <strong>n </strong>denotes the (unit) outward normal vector to Ω, and <em>f </em>satisfies the compatibility condition:

<em>.</em>

<ul>

 <li>Develop a second-order accurate FD method for solving this equation using the fictitious point method for deriving the approximations ∇<sup>2</sup><em>u </em>on the boundary (this will be exactly the same as we did in the 1-D problem). Implement a function similar to the fd2poissondst that solves the resulting system using the discrete cosine transform (DCT) and name the function fd2poissondct. Note that for this function it is unnecessary to pass in a value for <em>g </em>on the boundary. Also, the solution to this Poisson equation is unique up to an additive constant. As in the problem 5 of the last homework assignment, you will find that this property means the linear system from the FD2 approximation is singular with a vector of all ones (or non-zero scaled multiple) being the only eigenvector corresponding to zero eigenvalue. You need to deal with this singularity using the same approach as problem 5 of the previous homework assignment. You should set the arbitrary constant to zero.</li>

</ul>

Please download the updated dct and idct functions for solving this problem. Turn in a listing of your code.

2

<ul>

 <li>Use you code from part (a) to solve the Poisson equation with <em>f</em>(<em>x,y</em>) = −8<em>π</em><sup>2 </sup>cos(2<em>πx</em>)cos(2<em>πy</em>). An exact solution to the Poisson equation with this <em>f </em>is <em>u</em>(<em>x,y</em>) = cos(2<em>πx</em>)cos(2<em>πy</em>). Plot the difference between your computed solution and the exact solution for <em>m </em>= 2<sup>6 </sup>− 1. Also, generate a table showing the convergence of your solution to the true solution for <em>m </em>= 2<em><sup>k </sup></em>− 1, <em>k </em>= 4<em>,</em>5<em>,…,</em> The table should show the relative 2-norm of the error. Verify that your method is second-order accurate.</li>

</ul>

<ol start="4">

 <li><em>Implicit FD methods</em>

  <ul>

   <li>Using the technique from problem 4 of homework 2, derive the following implicit (compact) fourth-order accurate approximation to the 2-D Poisson equation <em>u<sub>xx </sub></em>+ <em>u<sub>yy </sub></em>= <em>f</em>:</li>

  </ul></li>

</ol>

<em>.</em>

<ul>

 <li>Write a code for solving the Poisson equation from problem 2, but now using the ninenode, compact, fourth-order accurate FD scheme. You can solve the problem in one of the following two ways. (1) use a sparse matrix library in whatever language you are using to implement the function to construct and solve the linear system as you did in part (a) of problem 2. (2) use SOR to solve the problem similar to what you did in problem 1 (again your code should be matrix-free). Use your code to solve the Poisson equation from problem 2 for various values of <em>m </em>and produce plots and tables that clearly show the fourth order accuracy of the method. Turn in a printed copy of your code and e-mail it to me.</li>

 <li>Extra credit (10 points): Do the same thing as 4b, but use the DST to solve the resulting</li>

</ul>

linear system. Produce plots showing your code gives the same results as the direct sparse solver. Illustrate the improved efficiency of the method by comparing the computational time it takes to compute a solution verses the direct sparse solver for various <em>m</em>.