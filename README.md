# Julia scripts

My Julia scripts used for numerical computation. I like [chaotic systems](https://en.wikipedia.org/wiki/Chaos_theory) and [Sturm-Liouville problems](https://en.wikipedia.org/wiki/Sturm%E2%80%93Liouville_theory) the most, so they are the main ones I solve in it.

## [Double pendulum](Double-pendulum.jl)

This script, as per its name, solves the problem of the double pendulum. Its differential function is largely copied from https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fau.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F27212%2Fversions%2F2%2Fcontents%2Fdouble_pendulum_ODE.m&embed=web (largely to save me a lot of time converting what equations I had from the Euler-Lagrange equations into a problem I could plug into ode78) therefore that segment of the script (which has been significantly modified to the point I don't think it even meets the threshold of originality to need this license) may fall under this license:

> Copyright (c) 2010, Alexander Erlich
> All rights reserved.

> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:

> * Redistributions of source code must retain the above copyright notice, this
>  list of conditions and the following disclaimer.

> * Redistributions in binary form must reproduce the above copyright notice,
>  this list of conditions and the following disclaimer in the documentation
>  and/or other materials provided with the distribution
>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
>AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
>IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
>DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
>FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
>DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
>SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
>CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
>OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
>OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

![Figure 1](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%201:%20theta2%20vs%20theta1.svg)
**Figure 1 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**

![Figure 2](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%202:%20dtheta1%20vs%20theta1.svg)
**Figure 2 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**

![Figure 3](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%203:%20dtheta2%20vs%20theta2.svg)
**Figure 3 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**

![Figure 4](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%204:%20dtheta2%20vs%20dtheta1.svg)
**Figure 4 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**

![Figure 5](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%205:%20theta1%20vs%20t.svg)
**Figure 5 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**

![Figure 6](https://fusion809.github.io/images/Graphs/Double-pendulum/Figure%206:%20theta2%20vs%20t.svg)
**Figure 6 (as of commit d43a31d97ec05ff896dc6246f742f8309c184a9e)**