/*
Copyright 2020 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code I've found for sampling from a Normal (and exponential) Distribution, and for the Normal CDF

///////////////////////////

/*  
 * Based on: M. M. Shepherd and J. G. Laframboise, "Chebyshev Approximation of 
 * (1+2x)exp(x^2)erfc x in 0 <= x < INF", Mathematics of Computation, Vol. 36,
 * No. 153, January 1981, pp. 249-253.
 *
 */  

float my_erfcxf (float x)
{
    float a, d, e, m, p, q, r, s, t;

    a = fmaxf (x, 0.0f - x); // NaN-preserving absolute value computation

    /* Compute q = (a-2)/(a+2) accurately. [0,INF) -> [-1,1] */
    m = a - 2.0f;
    p = a + 2.0f;
#if FAST_RCP_SSE
    r = fast_recipf_sse (p);
#else
    r = 1.0f / p;
#endif
    q = m * r;
    t = fmaf (q + 1.0f, -2.0f, a); 
    e = fmaf (q, -a, t); 
    q = fmaf (r, e, q); 

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
    p =              0x1.f10000p-15f;  //  5.92470169e-5
    p = fmaf (p, q,  0x1.521cc6p-13f); //  1.61224554e-4
    p = fmaf (p, q, -0x1.6b4ffep-12f); // -3.46481771e-4
    p = fmaf (p, q, -0x1.6e2a7cp-10f); // -1.39681227e-3
    p = fmaf (p, q,  0x1.3c1d7ep-10f); //  1.20588380e-3
    p = fmaf (p, q,  0x1.1cc236p-07f); //  8.69014394e-3
    p = fmaf (p, q, -0x1.069940p-07f); // -8.01387429e-3
    p = fmaf (p, q, -0x1.bc1b6cp-05f); // -5.42122945e-2
    p = fmaf (p, q,  0x1.4ff8acp-03f); //  1.64048523e-1
    p = fmaf (p, q, -0x1.54081ap-03f); // -1.66031078e-1
    p = fmaf (p, q, -0x1.7bf5cep-04f); // -9.27637145e-2
    p = fmaf (p, q,  0x1.1ba03ap-02f); //  2.76978403e-1

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    d = a + 0.5f;
#if FAST_RCP_SSE
    r = fast_recipf_sse (d);
#else
    r = 1.0f / d;
#endif
    r = r * 0.5f;
    q = fmaf (p, r, r); // q = (p+1)/(1+2*a)
    t = q + q;
    e = (p - q) + fmaf (t, -a, 1.0f); // residual: (p+1)-q*(1+2*a)
    r = fmaf (e, r, q);

    if (a > 0x1.fffffep127f) r = 0.0f; // 3.40282347e+38 // handle INF argument

    /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
    if (x < 0.0f) {
        s = x * x;
        d = fmaf (x, x, -s);
        e = expf (s);
        r = e - r;
        r = fmaf (e, d + d, r); 
        r = r + e;
        if (e > 0x1.fffffep127f) r = e; // 3.40282347e+38 // avoid creating NaN
    }
    return r;
}

///////////////////////////

