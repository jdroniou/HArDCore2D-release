# HArDCore2D
Hybrid Arbitrary Degree::Core 2D - Library to implement numerical schemes with edge and cell polynomial unknowns on generic 2D polygonal meshes.

The complete documentation is available here:

https://jdroniou.github.io/HArDCore2D-release/

See also the landing page for the HArDCore project: https://github.com/jdroniou/HArDCore/

The purpose of HArD::Core2D is to provide easy-to-use tools to code hybrid schemes, such as the Hybrid High-Order method. The data structure is described using intuitive classes that expose natural functions we use in the mathematical description of the scheme. For example, each mesh element is a member of the class 'Cell', that gives access to its diameter, the list of its edges (themselves members of the class 'Edge' that describe the geometrical features of the edge), etc. Functions are also provided to compute the key elements usually required to implement hybrid schemes, such as mass matrices of local basis functions, stiffness matrices, etc. The approach adopted is that of a compromise between readibility/usability and efficiency. 

As an example, when creating a mass matrix, the library requires the user to first compute the quadrature nodes and weights, then compute the basis functions at these nodes, and then assemble the mass matrix. This ensures a local control on the required degree of exactness of the quadrature rule, and also that basis functions are not evaluated several times at the same nodes (once computed and stored locally, the values at the quadrature nodes can be re-used several times). Each of these steps is however concentrated in one line of code, so the assembly of the mass matrix described above is actually done in three lines:

```
QuadratureRule quadT = generate_quadrature_rule(T, 2*m_K);<br>
boost::multi_array<double, 2> phiT_quadT = evaluate_quad<Function>::compute(basisT, quadT);<br>
Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, quadT);
```

Note that the `ElementQuad` class offers a convenient way to compute and store the quadrature rules and values of basis functions at the nodes, and makes it easy to pass these data to functions. More details and examples are provided in the documentation.

The implementations in this library follow general principles described in the appendix of the book "*The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications*" (D. A. Di Pietro and J. Droniou. 2019, 516p. url: https://hal.archives-ouvertes.fr/hal-02151813). High-order methods with hybrid unknowns have certain specificities which sometimes require fine choices, e.g. of basis functions (hierarchical, orthonormalised or not), etc. We refer to the aformentioned manuscript for discussion on these specificities. If using the code provided here, or part thereof, for a scientific publication, please refer to this book for details on the implementation choices.


This library was developed with the direct help and indirect advice of several people. Many thanks to them: Daniel Anderson, Lorenzo Botti, Daniele Di Pietro, Lachlan Grose, Tom Lemaitre, Liam Yemm.

The development of this library was partially supported by Australian Government through the Australian Research Council's Discovery Projects funding scheme (project number DP170100605).


