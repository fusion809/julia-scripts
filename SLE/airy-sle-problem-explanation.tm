<TeXmacs|1.99.8>

<style|generic>

<\body>
  <\hide-preamble>
    \;

    <assign|table-text|<macro|<localize|Table> 1: Selected eigenvalues>>
  </hide-preamble>

  <with|font-base-size|20|<\math>
    -<frac|<math-up|\<b-up-d\>><rsup|2>y<rsub|n><swell|>|<with|font-shape|italic|\<b-up-d\>><text|<with|font-series|medium|\<b-x\>>><rsup|2>>+\<b-x\>y<rsub|n>=\<lambda\><rsub|n>y<rsub|n>
  </math>>

  on <math|0\<leqslant\>x\<leqslant\>\<infty\>>, with
  <math|y<around*|(|0|)>=y<around*|(|\<infty\>|)>=0>.

  To this equation the solution is:

  <\math>
    y<rsub|n>=C<rsub|1 >Ai<around*|(|x-\<lambda\><rsub|n>|)>
  </math>

  Where <math|C<rsub|1>> is an arbitrary constant, and
  <math|-\<lambda\><rsub|n>> are the zeros of the Airy function,
  <math|Ai(x)>, on the negative x axis, i.e.
  <math|Ai<around*|(|-\<lambda\><rsub|n>|)>=0><space|1em><math|\<forall\>>
  <math|n>. <small-table*|<tabular|<tformat|<cwith|1|-1|1|-1|cell-background|>|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|1ln>|<cwith|1|-1|1|-1|cell-rborder|1ln>|<twith|table-halign|l>|<table|<row|<cell|<math|n>>|<cell|<math|\<lambda\><rsub|n>>>>|<row|<cell|1>|<cell|2.33811>>|<row|<cell|2>|<cell|4.08795>>|<row|<cell|3>|<cell|5.52056>>|<row|<cell|4>|<cell|6.78671>>|<row|<cell|5>|<cell|7.94413>>|<row|<cell|6>|<cell|10.0402>>>>>|>

  To numerically integrate this solution let us use a Chebyshev spectral
  method. If <math|t<rsub|m>> comprises points on the Chebyshev extrema grid,
  that is:

  <math|t<rsub|m>=-cos<around*|(|<frac|\<pi\>m|N>|)>>,<space|2em><math|m\<in\><around*|[|0,1,2,\<ldots\>,N|]>>

  then:

  <\math>
    x<rsub|m>=L<around*|(|<frac|1+t<rsub|m>|1-t<rsub|m>>|)>
  </math>

  .

  This is a rational transformation from the Chebyshev extrema grid to the
  semi-infinite domain.\ 
</body>

<\initial>
  <\collection>
    <associate|font-base-size|20>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|?|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal||<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>