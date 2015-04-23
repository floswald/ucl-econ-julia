Finite Horizon Life Cycle model with AR(1) iid income process (agent chooses C and A1, and might
retire at a given age exogenously)
------------------------------------------------------------------------------------------------------

This code gets the same policy rules as Cormac & Monica Costa version in matlab (Dynamic Economics
28-29 October 2013, IFS). There are some minor changes on the algorithm, but is as similar as possible
to make it clear how to translate models from Matlab into Julia. There are many functions that were
not translated as there are native versions of it on the basic packages.

For sure there are still bugs or even untranslated functions in some switchs;
so just let me know if you find them. Folder "matlabObj" has exported output from matlab in order to
compare using the exact policyA1 and shocks that are generated there.

Just include main_v5 and it will load all the required files. I've ordered in folders as it makes it
more organized
Graphics are exported into "output/images", using Gadfly and Cairo. As it is slow, it is within an
impossible condition (1==0) that you need to deactivate if you want to generate those graphs.
Alternatively, you can use the notebook for IJulia (Graphs Main5.ipynb) with some interactive features

Paul

