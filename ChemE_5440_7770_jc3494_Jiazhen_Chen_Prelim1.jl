### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 9811211c-f893-40ae-bf79-b4a6c274d3e1
# Load essential packages
begin

	# external packages -
	using PlutoUI
	using Plots
	using PrettyTables
	using BlackBoxOptim
end

# ╔═╡ 7b5a2550-aae9-11ec-1d1a-5978db55a3c7
md"""
#### **Prelim 1**
##### **Net ID:** jc3494 **Student Name:** Jiazhen Chen
"""

# ╔═╡ c64a141f-6610-4edd-86d0-1d8056213de4
md"""
##### **Question 1:**
##### **Part a)**
From the question prompt, we have the following equations:

$$\begin{eqnarray}
𝓖_{j}+𝓡_{X} &\rightleftharpoons& {(𝓖_{j}:𝓡_{X})_{O}}\\
{(𝓖_{j}:𝓡_{X})_{C}} &\longrightarrow& {(𝓖_{j}:𝓡_{X})_{O}}\\
{(𝓖_{j}:𝓡_{X})_{O}} &\longrightarrow& {𝓖_{j}+𝓡_{X}+𝓶_j^*}\\
{(𝓖_{j}:𝓡_{X})_{O}} &\longrightarrow& {𝓖_{j}+𝓡_{X}+𝓶_j}\\
𝓡_{X,j} & = & 𝓚_{E,j}(𝓖_{j}:𝓡_{X})_{O}
\end{eqnarray}$$ 

Assume the steady state assumption for the intermediates:


$$\frac{d(𝓖_{j}:𝓡_{X})_{C}}{dt} = k_{+j}𝓖_{j}𝓡_{X} - k_{-j}(𝓖_{j}:𝓡_{X})_{C}\simeq{0}$$

and

$$\frac{d(𝓖_{j}:𝓡_{X})_{O}}{dt} = k_{I,j}(𝓖_{j}:𝓡_{X})_{C} - k_{A,j}(𝓖_{j}:𝓡_{X})_{O} - k_{E,j}(𝓖_{j}:𝓡_{X})_{O}\simeq{0}$$

From the above equations for steady-state assumption, we can rewrite the concentration for the intermediates using $𝓖_{j}$, $𝓡_{X}$, and the rate constants. 

$(𝓖_{j}:𝓡_{X})_{C} = \frac{k_{+j}}{k_{-j}+k_{I}}𝓖_{j}𝓡_{X}$
$(𝓖_{j}:𝓡_{X})_{O} = \frac{k_{I,j}}{k_{A,j}+k_{E,j}}(𝓖_{j}:𝓡_{X})_{C}$

The two equations above can then be further reduced by using the definitions of saturation and time constants.

$(𝓖_{j}:𝓡_{X})_{C} = \frac{1}{k_{X,j}}𝓖_{j}𝓡_{X}$
$(𝓖_{j}:𝓡_{X})_{O} = \frac{1}{k_{X,j} τ_{X,j}}𝓖_{j}𝓡_{X}$

Similarly, for substrates $i≠j$

$(𝓖_{i}:𝓡_{X})_{C} = \frac{1}{k_{X,i}}𝓖_{i}𝓡_{X}$
$(𝓖_{i}:𝓡_{X})_{O} = \frac{1}{k_{X,i} τ_{X,i}}𝓖_{i}𝓡_{X}$

With the relations, the rate equation for gene $j$ can be rewritten as:

$𝓻_{X,j}=k_{E,j}\frac{1}{k_{X,j} τ_{X,j}}𝓖_{j}𝓡_{X}$

We can solve for $𝓡_{X}$ by substituting the four equations above into the equation for the balance of $RNAP$.

$𝓡_{X,T} = 𝓡_{X} + (𝓖_{j}: 𝓡_{X})_{C} + (𝓖_{j}: 𝓡_{X})_{O} + \sum_{i=1,j}^{\mathcal{N}}\{(𝓖_{i}: 𝓡_{X})_{C}+(𝓖_{j}: 𝓡_{X})_{O}\}$

Solving for $𝓡_{X}$:

$𝓡_{X} = \frac{𝓡_{X,T}}{1+k_{X,j}^{-1}𝓖_{j}+k_{X,j}^{-1}τ_{X,j}^{-1}𝓖_{j}+\sum_{i=1,j}^{\mathcal{N}}\frac{1+τ_{X,i}}{k_{X,i}τ_{X,i}}𝓖_{i}}$

Replacing $𝓡_{X}$ in the rate equation with the *RHS* above we arrive at:

$𝓻_{X,j}=k_{E,j}𝓡_{X,T}(\frac{𝓖_{j}}{k_{X,i} τ_{X,i}+(1+ τ_{X,i})𝓖_{j}+𝓔_{j}})$

where

$𝓔_{j} = \sum_{i=1,j}^{\mathcal{N}}\frac{k_{X,j}τ_{X,j}}{k_{X,i}τ_{X,i}}(1+τ_{X,i})𝓖_{i}$ 

Thus the equation is proven.\


##### **Part b)**
As defined by the problem, when neglecting RNAP abort rate we have:

$τ_{X,j} = \frac{k_{A,j}+k_{E,j}}{k_{X,j}}\simeq\frac{k_{E,j}}{k_{X,j}}$

When $τ_{X,j}<<1$,  we can conclude that $k_{E,j}<<k_{X,j}$, and therefore transcription is **elongation limited**. 
\
\
\
"""

# ╔═╡ 66f6bf9a-a606-4a13-9313-7e4e63934839
md"""
##### **Question 2:**
##### **Part a)**

The three different micro-states for the model is defined as follows:
* __State $s = 0$__: **No substrate** and **no effector** is bound, **PFK** is at base state. No reaction is possible
* __State $s = 1$__: **Substrate** is bound to **PFK**, but **no effector** is bound. The reaction is possible but the rate is slow.
* __State $s = 2$__: The **substrates and effector** are all bound to **PFK**. The reaction is occurring at a high rate.

Assume that for state **$s = 0$** and **$s = 1$**, $f_0 = f_1 = 1$.


For state **$s = 2$**, $f_3$ satisfies the hill type binding function:

$$f_{2} = \frac{(AMP/K_{2})^{n_{2}}}{(1+(AMP/K_{2})^{n_{2}})}$$

Following the question prompt, we have the following equations:

$$p_{i} = \frac{W_{i}f_{i}}{Z}$$
where

$Z=\sum_{i=0}^{\mathcal{2}}W_{i}f_{i}$

States $s = 1$ and $s = 2$ are microstates that are able to catalyze the reaction, thus the control function $v$ can be written as:

$$v = \frac{W_{1}f_{1}+W_{2}f_{2}}{Z}$$

The reaction rate can be written as:

$\hat{r}_{1} = r_{1}v$

where

$r_{1} = k_{cat}E_{1}(\frac{F6P}{K_{F6P}+F6P})(\frac{ATP}{K_{ATP}+ATP})$
"""

# ╔═╡ 27cd2930-de3e-486c-aeeb-3f3fe2377a56
md"""
##### **Part b)**

The above setup is defined as a function and all unknown variables are estimated by using the *BlackBoxOptim.jl* optimizer package, and the experimental data provided in **Table 1**.

"""

# ╔═╡ 2086a2af-3d4b-42ae-9b9b-31509d4e4541
begin
	#This section converts all known values and all experimental values into vector form.
	
	#input all experimental values
	AMP = [0.0 0.055 0.093 0.181 0.405 0.990 1.0] #unit mM
	rhat = [3.0 6.3 29.8 52.0 60.3 68.7 68.9] #unit uM/h
	std = [0.59 1.20 5.7 10.2 11.8 13.3 10.0] #uniut uM/h
	
	#calculate the maximum rate
	kcat = 0.4 #unit s-1
	E1 = 0.12 #unit uM
	F6P = 0.1 #unit mM
	ATP = 2.3 #unit mM
	k_F6P = 0.11 #unit mM
	k_ATP = 0.42 #unit mM
	r1 = kcat*E1*(F6P/(F6P+k_F6P))*(ATP/(ATP+k_ATP))*3600 #unit uM/h
	
	nothing
end

# ╔═╡ abf220a8-5146-4e90-8f73-cdee3f51136d
md"""

I have decided to use the **_BlackBoxOptim.jl_** because, unlike **_Optim.jl_**,  it contains algorithms for global optimization, suitable for poorly behaving problems with multiple local minimums, such as this one.  

The error function that I have tried to minimize is based on a **"weighted least squares" method**. Using this method, the seven experimental data points contributes differently depending on their standard deviation:

$Error=\sum_{i=1}^{\mathcal{7}}\frac{E_{i}^2}{σ_{i}^2}$

where $E_{i}$ is the difference between modeled and experimental results.

The input variables are the unknowns in the model and are held in the vector $\vec{v}$. For simplicity, **$W_{0}$ was set to 1** by definition to represent the accessibility of the base state. 

The optimization range was selected based on the following principles, and arbitrary values were selected: 

* All parameters should be within $[0, ∞)$. 
*  $K$ is typically smaller than 1 mM (checked on BioNumbers).
*  $n_{2}$ should be $O(1)$.  
*  $W$ could be very large.

The code used is shown below:

"""

# ╔═╡ fe9951d0-feeb-4308-aa73-879771bfebdc
begin
	function f(v) 
		
	#values of fi
    x = (AMP ./v[1]) .^v[2]

	f_array = ones(3,7) 

    f_array[3,:] = x./(1 .+x) 
    	
	#values of W
	W_array = [
		1
		v[3]
		v[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1 ./Z).*microstate_array
	v = p_array[2,:] .+p_array[3,:] 

	#compute error function
	err= r1.*v .- vec(rhat) 
	w = 1 ./(vec(std)).^2 
	abserr = err.^2 .*w
	
	return sum(abserr)
	end

	nothing
end
	

# ╔═╡ 1e42f1f9-234d-4412-856d-235b07348882
# This part is for setting up the optimization of parameters for model fitting. A max time of 15s is typically sufficient for search convergence (stop due to function evaluation), and obtaining reproducible results.

begin
	res = bboptimize(f; SearchRange = [(0.0, 5),(0.0, 10.0),(0,1e9),(0,1e9)], Method = :adaptive_de_rand_1_bin, MaxNumStepsWithoutFuncEvals = 100, MaxTime = 15.0)
	vars = best_candidate(res)
	nothing
end

# ╔═╡ f77e448d-31a9-46f2-ab8e-b1cb66249968
# Show optimizer result. stop_reason should be "Too many steps(101) without...". In rare cases it is not, please rerun the previous cell.
res

# ╔═╡ 6bac3eb6-8a92-4acd-9642-6acf8d27610e
md"""
The calculated and predetermined coefficients are summarized in the table below. 
"""

# ╔═╡ 73029713-0b7d-4f9d-a934-42385e2b6afc
with_terminal() do

	# initialize -
	state_table = Array{Any,2}(undef,5,2)
	
	par_name = ["K₂","n₂","W₀","W₁","W₂"]
	par_val = [vars[1],vars[2],1, vars[3], vars[4]]
	
	# populate the state array -
	for n ∈ 1:5

		state_table[n,1] = par_name[n]
		state_table[n,2] = par_val[n]

	end

	# header -
	header_row = (["Parameters","Values"])
	pretty_table(state_table;header=header_row,title = "Best Fit Parameters for Q2")
	
end

# ╔═╡ 1c73aca2-2862-4979-9bcc-067ad8e00836
md"""
##### **Part c)**
"""

# ╔═╡ 568ad739-417d-4492-9122-a4b22f042db0
md"""
The parameters are substituted back into the model, and the resulting reaction rate for any *AMP* concentration could be determined.
"""

# ╔═╡ e2e342e5-9839-48a0-8d1e-49cac2c237da
#Function for finding the rate using optimized parameters.
begin
	function rate(x)
	f_array = [
		1
		1
		((x/vars[1])^vars[2])/(1+(x/vars[1])^vars[2])
	]

	
	#values of W
	W_array = [
		1
		vars[3]
		vars[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1/Z).*microstate_array
	v = p_array[2]+p_array[3]

	#calculate rate
	model_rate = r1.*v
	
	return model_rate
	end
	
	nothing
end

# ╔═╡ 3f3bef97-c340-47b6-b02d-ddd703a1c09f
#altered error function, this section will be used for plotting the dashed line in grey.

begin
	
	function f2(v) 
		
	#values of fi
    x = (AMP ./v[1]) .^v[2]

	f_array = ones(3,7) 

    f_array[3,:] = x./(1 .+x) 
    	
	#values of W
	W_array = [
		1
		v[3]
		v[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1 ./Z).*microstate_array
	v = p_array[2,:] .+p_array[3,:] 

	#compute error function
	err= r1.*v .- vec(rhat) 
	w = 1 ./(vec(std)).^2 
	abserr = err.^2
	
	return sum(abserr)
	end

	res2 = bboptimize(f2; SearchRange = [(0.0, 5),(0.0, 10.0),(0,1e9),(0,1e9)], Method = :adaptive_de_rand_1_bin, MaxNumStepsWithoutFuncEvals = 100, MaxTime = 15.0)
	vars2 = best_candidate(res2)
	
	function rate2(x)
	f_array = [
		1
		1
		((x/vars2[1])^vars2[2])/(1+(x/vars2[1])^vars2[2])
	]

	
	#values of W
	W_array = [
		1
		vars2[3]
		vars2[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1/Z).*microstate_array
	v = p_array[2]+p_array[3]

	#calculate rate
	model_rate = r1.*v
	
	return model_rate
	end
	nothing	
end

# ╔═╡ 8a208af2-38c1-402f-bda6-a954089ddc72
begin

	# get MM parameter values -

	number_of_steps = 100

	# substrate range -
	AMP_array = range(0.0,stop=1.0,length=number_of_steps) |> collect;

	# initialize space -
	v_array = Array{Float64,1}(undef,number_of_steps)
	v_array2 = Array{Float64,1}(undef,number_of_steps)

	#compute the rate -
	for i ∈ 1:number_of_steps
		
		v_array[i] = rate(AMP_array[i])
		v_array2[i] = rate2(AMP_array[i])
		
	end

	plot(AMP_array,v_array,xlims=(0.0,1.0),ylims=(0.0,90), lw = 2, linecolor = "black", label = "Best Fit Weighted Sum Model")
	plot!(AMP_array,v_array2,xlims=(0.0,1.0),ylims=(0.0,90), lw = 1, linecolor = "grey", linestyle = :dash, label = "Best Fit Sum Model")
	
	scatter!(AMP,rhat,label = ["Experimental Value" "" "" "" "" "" ""], marker = (:diamond, :black),yerror = std)
	xlabel!("AMP Concentration (mM)",fontsize=18)
	ylabel!("Rate v (μM/h)",fontsize=18)
	plot!(plot_title = "Comparison of Best Fit Model with Experimental Results", plot_titlefontsize = 12)
	
end

# ╔═╡ df86bb0a-fd21-4a96-8312-54e015aac475
md"""
The **black solid line** shows the best fit model according to my **weighted sum** error function. 

As can be seen, the model had been effective in predicting the general trend but failed to capture the specific shape of the curve. The "S" shape of the experimental results was more subtle than what the model showed. In addition, when the model is best-fit, it also was not able to match well with the assumption that at infinite *AMP* concentration, $v=v_{max}$.

Better fitting for the latter points could be easily achieved by changing the error function, for example, using **_absolute sum (grey dashed line)_** (cell containing code hidden directly above the figure). 

However, doing so would result in the model being inaccurate for other points, in this case, for low concentrations. This seems to be an issue with the model itself and could not be resolvable by parameter tuning.
\
\
\
"""

# ╔═╡ 981e54a3-9198-4ff5-8af6-37e4e6ca5718
md"""
##### **Question 3:**
##### **Part a)**
Conversion is done using the equation:

$\frac{C_{mRNA}\times V_{c}\times 10^{-15}}{M_{c}}$

where $C_{mRNA}$ is in $[nmol/L]$, $V_{c}$ is in $[μm^3]$, and $M_{c}$ is in $[gDW]$.
"""

# ╔═╡ 6a0d60e8-423d-4a21-a9f1-553b16d1b59b
#Convert Question Prompt to Code
begin
	conc_mRNA = [19, 21, 41, 67, 86, 93, 93] # nmol/L
	Vc = 2.75 # um^3
	Mc = 2e-13 # gDW
	fac = 1e-15 # conversion factor between L and um^3
	IPTG = [0.0 5e-4 0.005 0.012 0.053 0.216 1.0] #unit mM
	
	#calculate nmol/gDW
	n_array = (fac*conc_mRNA*Vc)./Mc
	
	nothing
end

# ╔═╡ 73a56522-388f-41a9-aa48-f88f996ed0e9
md"""
Calculated values are summarized in the table below.
"""

# ╔═╡ 24ba0aa1-3d08-49bc-9e1d-8b104931c9c9
with_terminal() do

	# initialize -
	n_table = Array{Any,2}(undef,7,2)
	
	
	# populate the state array -
	for n ∈ 1:7

		n_table[n,1] = IPTG[n]
		n_table[n,2] = n_array[n]

	end

	# header -
	header_row = (["IPTG","<n> nmol/gDW"])
	pretty_table(n_table;header=header_row,title = "Converted <n> Values")
	
end

# ╔═╡ a07eb6b6-bc71-488a-83e2-c43acf8a5fa4
md"""
##### **Part b)**
For the steady-state approximation, we can assume that the change of mRNA can be neglected, resulting in the steady-state mRNA balance:

$\dot{m_{i}}=r_{X,i}\bar{u}_{i}-(μ+θ_{m,i})m_{i}\simeq0$

By simplification of the above equation:

$m^*=\frac{r_{X,i}\bar{u}_{i}}{μ+θ_{m,i}}=K_{X}⋅\bar{u}_{i}$

In this case, $K_{X}$ can be understood as the max mRNA abundance and can be estimated from the final values of the table. $K_{X}$ has the units nM/gDW, and is thus equal to the largest mass-specific concentration:

$K_{X}≡1.279$

The model setup of $\bar{u}$ in  Question 3 is analogous to the control function $v$ as in Question 2. The data shows that with no *IPTG* induction, transcription was still possible. Thus:

* __State $s = 0$__: Gene $G_{i}$ is at base state. No reaction is possible
* __State $s = 1$__: *RNAP* is bound to $G_{i}$, but no *IPTG* is available. Transcription possible but very slow.
* __State $s = 2$__: *RNAP* and *IPTG* is bound to $G_{i}$. Transcription is occurring at high rate.

Following the exact logic as in Question 2:

Assume that for state **$s = 0$** and **$s = 1$**, $f_0 = f_1 = 1$.


For state **$s = 2$**, $f_3$ satisfies the hill type binding function:

$$f_{2} = \frac{(I/K_{2})^{n_{2}}}{(1+(I/K_{2})^{n_{2}})}$$

Following the question prompt, we have the following equations:

$$p_{i} = \frac{W_{i}f_{i}}{Z}$$
where

$Z=\sum_{i=0}^{\mathcal{2}}W_{i}f_{i}$

and

$W_{0} ≝ 1$

States $s = 1$ and $s = 2$ are microstates that are able to transcribe, thus the promoter activity function can be written as:

$$\bar{u}_{i} = \frac{W_{1}f_{1}+W_{2}f_{2}}{Z}$$

"""

# ╔═╡ b2ae2b86-6b1c-45d6-882a-6e8fb2e5edbf
md"""
##### **Part c)**
The growth and decay rates are calculated using:

$μ = \frac{log(2)}{τ_{d}}$

and

$θ_{m,i} = \frac{log(2)}{τ_{h}}$

The limit for the rate of transcription $r_{X,i}$ is calculated by:

$r_{X,i} = V_{max,i}(\frac{G_{i}}{K_{X,i}+G_{i}}) = K_{X}(μ+θ_{m,i})$

with

$V_{max,i}=k_{3,i}\left[RNAP\right]_{T}$
$k_{3,i} = \frac{e_{X}}{L_{lacZ}}$

Additionally, $K_{X,i}$ can be solved from the above equations. Using the definition of $K_{X,i}$ and assuming $k_{D,i}=550 nM$, other kinetic parameters can be calculated, although it is unnecessary for this problem:

$K_{X,i}=\frac{k_{2,i}+k_{3,i}}{k_{1,i}}$
$k_{D,i}= k_{2,i}/k_{1,i}$
 
"""

# ╔═╡ 7dc0d32f-e164-4c2c-94a1-5d1c326d4a56
#Converting all data to code:
begin
	Tdouble = 40*60 #cell doubling time; unit sec
	Tdeg = 5*60 #mRNA degradation time; unit sec
	ex = 35 #extension rate; unit nt/s
	L = 3075 #lacZ length; unit nt
	G = 2 # unit copies/cell or nM
	RNAP = 4600 # unit copies/cell or nM
	RNAP_DW = (fac*RNAP*Vc)./Mc #unit nM/gDW (converted to dry weight unit for converted comparison)
	
	
	mu = log(2)/Tdouble #unit 1/sec
	theta_mi = log(2)/Tdeg # unit 1/sec
	K_D = 550 # unit nM BNID:103590 
	k_3 = ex/L # unit 1/sec
	
	Vmaxi = k_3*RNAP_DW # unit nM/(gDW⋅sec)
	K_X =  1.279 #unit nM/gDW
	
	#calculating range of dataset
	H = [20 26 44 69 88 95 94]
	L = [18 17 37 65 84 91 92]

	H_DW = vec((fac*H*Vc)./Mc)
	L_DW = vec((fac*L*Vc)./Mc)
	
	R = H_DW .-L_DW
	
	nothing
end

# ╔═╡ 2a25a071-e439-42ba-a278-785ef64c363f
begin
	function g(v) 
		
	#values of fi
    x = (IPTG ./v[1]) .^v[2]

	f_array = ones(3,7) 

    f_array[3,:] = x./(1 .+x) 
    	
	#values of W
	W_array = [
		1
		v[3]
		v[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1 ./Z).*microstate_array
	u_tilde = p_array[2,:] .+p_array[3,:] 

	#compute error function
	err= K_X.*u_tilde .- vec(n_array) 
	w = 1 ./(vec(R)).^2 
	abserr = err.^2 .*w
	return sum(abserr)
	end
	nothing
end

# ╔═╡ 59c32f10-84cd-4af1-9b39-16e750efbbc9
md"""
The range of $K_{2}$ for *IPTG* and *RNAP* interaction cannot be found on Bionumbers. Therefore, the value is determined through trial and error, by adjusting until the optimal solution is bounded by the search range. Values were chosen arbitrarily. Design guidelines remained the same as **Question 2**, including:
* All parameters should be within $[0,∞)$.
*  $n_{2}$ should be of $O(1)$.
*  $W$ may be very large.
Furthermore, when defining the function, the weights are defined by the range(difference between upper and lower limit) as opposed to standard deviation in Question 2:

$Error=\sum_{i=1}^{\mathcal{7}}\frac{E_{i}^2}{Range_{i}^2}$

where Range is the difference between converted $high\ (mRNA/cell)$ and $low\ ( mRNA/cell)$.
"""

# ╔═╡ 3f5b279b-e2d4-48d1-ac2f-afc24b2c01b7
begin
	res3 = bboptimize(g; SearchRange = [(0.0,5.0),(0.0, 10.0),(0,1e9),(0,1e9)], Method = :adaptive_de_rand_1_bin, MaxNumStepsWithoutFuncEvals = 100, MaxTime = 20.0)
	vars3 = best_candidate(res3)
	nothing
end

# ╔═╡ 1a0cc0e9-d459-46be-8a80-0b4da76eb06e
# A MaxTime setting of 10s is typically enough for full convergence. If stop_reason does not return "Too many steps (101) without...", please run the previous cell again.
res3

# ╔═╡ 60a14e03-e988-4864-a655-ba7f26592583
with_terminal() do

	# initialize -
	state_table = Array{Any,2}(undef,5,2)
	
	par_name = ["K₂","n₂","W₀","W₁","W₂"]
	par_val = [vars3[1],vars3[2],1, vars3[3], vars3[4]]
	
	# populate the state array -
	for n ∈ 1:5

		state_table[n,1] = par_name[n]
		state_table[n,2] = par_val[n]

	end

	# header -
	header_row = (["Parameters","Values"])
	pretty_table(state_table;header=header_row,title = "Best Fit Parameters for Q3")
	
end

# ╔═╡ caa7e6ff-b7c6-4820-9afa-a80bb89f093e
md"""
##### **Part D)**
"""

# ╔═╡ 5b8cde28-c3c2-4398-8527-4dd3714c76a5
begin
function rate3(x)
	f_array = [
		1
		1
		((x/vars3[1])^vars3[2])/(1+(x/vars3[1])^vars3[2])
	]

	
	#values of W
	W_array = [
		1
		vars3[3]
		vars3[4]
	]

	#calculate factor v 
	microstate_array = f_array.*W_array
	Z = sum(microstate_array, dims=1)
	p_array = (1/Z).*microstate_array
	u_tilde = p_array[2]+p_array[3]

	#calculate rate
	model_rate = K_X.*u_tilde
	
	return model_rate
	end
	
	nothing
end

# ╔═╡ 44ae8fed-5700-4bfb-99bf-02b5a6720c30
begin
	
	# substrate range -
	IPTG_array = range(0.0,stop=1.0,length=number_of_steps) |> collect;

	# initialize space -
	vq3_array = Array{Float64,1}(undef,number_of_steps)
	
	#compute the rate -
	for i ∈ 1:number_of_steps
		
		vq3_array[i] = rate3(IPTG_array[i])
		
	end

	plot(IPTG_array,vq3_array,xlims=(0.0,1.0),ylims=(0.0,1.7), lw = 2, linecolor = "black", label = "Best Fit Weighted Sum", title = "Best Fit mRNA Concentration for Full Range IPTG",titlefontsize = 14)
	scatter!(vec(IPTG),n_array,label = ["Experimental Value" "" "" "" "" "" ""], 
	yerror = (n_array.-L_DW,H_DW.-n_array), marker = (:diamond, :black))
	xlabel!("IPTG Concentration (mM)",labelfontsize=10)
	ylabel!("mRNA Concentration (nM/gDW)",labelfontsize=10)


end

# ╔═╡ f405d547-488d-4471-9272-5f07d46df0ff
md"""
The best-fit modeled curve was plotted along with the experimental data. Error bars represent the high and low mRNA values. 

In general, the plot had a good fit and captured the general trend very well. At IPTG = 0.012 mM, the predicted value was outside the error range, as can be seen below.  
"""

# ╔═╡ bca58dbe-3bb7-4020-b5b7-f50a177d2e2a
begin
	plot(IPTG_array,vq3_array,xlims=(0.0,0.02),ylims=(0.0,1.7), lw = 2, linecolor = "black", label = "Best Fit Weighted Sum", title = "Best Fit mRNA Concentration at Low IPTG",titlefontsize = 14)
	scatter!(vec(IPTG),n_array,label = ["Experimental Value" "" "" "" "" "" ""], 
	yerror = (n_array.-L_DW,H_DW.-n_array), marker = (:diamond, :black))
	xlabel!("IPTG Concentration (mM)",labelfontsize=10)
	ylabel!("mRNA Concentration (nM/gDW)",labelfontsize=10)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BlackBoxOptim = "a134a8b2-14d6-55f6-9291-3336d3ab0209"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[compat]
BlackBoxOptim = "~0.6.1"
Plots = "~1.27.2"
PlutoUI = "~0.7.37"
PrettyTables = "~1.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BlackBoxOptim]]
deps = ["CPUTime", "Compat", "Distributed", "Distributions", "HTTP", "JSON", "LinearAlgebra", "Printf", "Random", "SpatialIndexing", "StatsBase"]
git-tree-sha1 = "41e347c63757dde7d22b2665b4efe835571983d4"
uuid = "a134a8b2-14d6-55f6-9291-3336d3ab0209"
version = "0.6.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CPUTime]]
git-tree-sha1 = "2dcc50ea6a0a1ef6440d6eecd0fe3813e5671f45"
uuid = "a9c8d775-2e2e-55fc-8582-045d282d599e"
version = "1.0.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "c43e992f186abaf9965cc45e372f4693b7754b22"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.52"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "90021b03a38f1ae9dbd7bf4dc5e3dcb7676d302c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpatialIndexing]]
git-tree-sha1 = "fb7041e6bd266266fa7cdeb80427579e55275e4f"
uuid = "d4ead438-fe20-5cc5-a293-4fd39a41b74c"
version = "0.1.3"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6976fab022fea2ffea3d945159317556e5dad87c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─7b5a2550-aae9-11ec-1d1a-5978db55a3c7
# ╟─c64a141f-6610-4edd-86d0-1d8056213de4
# ╟─66f6bf9a-a606-4a13-9313-7e4e63934839
# ╟─27cd2930-de3e-486c-aeeb-3f3fe2377a56
# ╠═2086a2af-3d4b-42ae-9b9b-31509d4e4541
# ╟─abf220a8-5146-4e90-8f73-cdee3f51136d
# ╠═fe9951d0-feeb-4308-aa73-879771bfebdc
# ╠═1e42f1f9-234d-4412-856d-235b07348882
# ╠═f77e448d-31a9-46f2-ab8e-b1cb66249968
# ╟─6bac3eb6-8a92-4acd-9642-6acf8d27610e
# ╟─73029713-0b7d-4f9d-a934-42385e2b6afc
# ╟─1c73aca2-2862-4979-9bcc-067ad8e00836
# ╠═568ad739-417d-4492-9122-a4b22f042db0
# ╠═e2e342e5-9839-48a0-8d1e-49cac2c237da
# ╟─3f3bef97-c340-47b6-b02d-ddd703a1c09f
# ╟─8a208af2-38c1-402f-bda6-a954089ddc72
# ╟─df86bb0a-fd21-4a96-8312-54e015aac475
# ╟─981e54a3-9198-4ff5-8af6-37e4e6ca5718
# ╠═6a0d60e8-423d-4a21-a9f1-553b16d1b59b
# ╟─73a56522-388f-41a9-aa48-f88f996ed0e9
# ╟─24ba0aa1-3d08-49bc-9e1d-8b104931c9c9
# ╟─a07eb6b6-bc71-488a-83e2-c43acf8a5fa4
# ╟─b2ae2b86-6b1c-45d6-882a-6e8fb2e5edbf
# ╠═7dc0d32f-e164-4c2c-94a1-5d1c326d4a56
# ╠═2a25a071-e439-42ba-a278-785ef64c363f
# ╠═59c32f10-84cd-4af1-9b39-16e750efbbc9
# ╠═3f5b279b-e2d4-48d1-ac2f-afc24b2c01b7
# ╠═1a0cc0e9-d459-46be-8a80-0b4da76eb06e
# ╟─60a14e03-e988-4864-a655-ba7f26592583
# ╟─caa7e6ff-b7c6-4820-9afa-a80bb89f093e
# ╟─5b8cde28-c3c2-4398-8527-4dd3714c76a5
# ╟─44ae8fed-5700-4bfb-99bf-02b5a6720c30
# ╟─f405d547-488d-4471-9272-5f07d46df0ff
# ╟─bca58dbe-3bb7-4020-b5b7-f50a177d2e2a
# ╠═9811211c-f893-40ae-bf79-b4a6c274d3e1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
