using JuMP


pi = 30.*rand(24,6)

ssl = 13;
hh = 24;
num_types = 15;

const pick_o = vec(ones(24));
pick_o[13:17] = 2;
pick_o[18:24] = 3;
pick_o[1:5] = 4;
const pick_fine = vec(ones(24));                                                                        
pick_fine[14:18]=2;                                                                                
const lease=vec(ones(24).*105);                                                                                
lease[5:16]=116;                                                                                    
const work_time = vec(ones(24));
work_time[5:17] = 2; 


breakm = rand((24,12,15,6));
active_stop = round(rand((24,12,15,6)));
active_nonstop = round(3.*rand((24,12,15,6)));
inactive_start = round(10.*rand((24,15)));
inactive_nonstart = round(10.*rand((24,15)));
pw = [0.00255578, 0.0886157, 0.408828, 0.408828, 0.0886157, 0.00255578];
nwage = length(pw);
constr_mat = round(rand((num_types,24)))

m = Model()

@defVar(m,  0.0001 <= p_est[1:hh,1:ssl-1,1:num_types,1:nwage] <= 0.9999);                             
@defVar(m,  0.0001 <= q_est[1:hh,1:num_types] <= 0.9999);
@defVar(m, -1000 <= EV[1:hh,1:ssl-1,1:num_types,1:nwage] <= 1700);
@defVar(m, -1000 <= EVA[1:hh,1:ssl-1,1:num_types] <= 1700);
@defVar(m,0 <= λ_1 <= 20);                                                                    
@defVar(m,0 <= λ_2 <= 20);                                                                    
@defVar(m,0 <= λ_0[1:length(unique(pick_o))] <= 250);                                                   
@defVar(m,-100 <= f[1:2] <= 400);
@defVar(m, 0.1 <= σ <= 80);                                                                      
@defVar(m, 0.1 <= σ_s <= 80);                                                               
@defVar(m,-500 <= μ[1:2] <= 600);                             

setValue(λ_1, -8.0);
setValue(λ_2, .8);
setValue(μ[1], 90.0);
setValue(μ[2], 90.0);


setValue(σ_s, 20.4);
setValue(σ, 20.0);


for h in [unique(pick_o)]
setValue(λ_0[h], 20.0)
end;

for h in [unique(pick_fine)]
    setValue(f[h], 30.0)  
end;

for h=1:hh,sl=1:ssl-1, c=1:num_types, a=1:nwage                 
    setValue(p_est[h,sl,c,a], 0.06*sl)  
end;

for h=1:hh,sl=1:ssl-1, c=1:num_types, a=1:nwage                             
setValue(EV[h,sl,c,a], 15)
end;

for h=1:hh,sl=1:ssl-1, c=1:num_types                                      
setValue(EVA[h,sl,c], 15)
end;

for h=1:hh, c=1:num_types                                                
    setValue(q_est[h,c], 0.5)hhend;


@setNLObjective(m, Max, sum{(1/10000)*(active_stop[h,sl,c,a]*log(p_est[h,sl,c,a])+active_nonstop[h,sl,c,a]*log(1-p_est[h,sl,c,a])),h=1:hh,sl=1:ssl-1,c=1:num_types, a=1:nwage}+sum{(1/10000)*(inactive_start[h,c]*log(q_est[h,c])+inactive_nonstart[h,c]*log(1-q_est[h,c])),h=1:hh,c=1:num_types})

for h=1:hh-1,sl=1:ssl-2, c=1:num_types, a=1:nwage
    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[h+1,sl+1,c])/σ)))
end


for h=1:hh-1,sl=ssl-1, c=1:num_types, a=1:nwage
    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
end


for h=hh,sl=1:ssl-2, c=1:num_types, a=1:nwage
    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[1,sl+1,c])/σ)))
end


for h=hh,sl=ssl-1, c=1:num_types, a=1:nwage
    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
end

for h=1:hh-1,sl=1:ssl-2, c=1:num_types, a=1:nwage
    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[h+1,sl+1,c])/σ)))
end

for h=1:hh-1,sl=ssl-1, c=1:num_types, a=1:nwage
    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
end

for h=hh,sl=1:ssl-2, c=1:num_types, a=1:nwage
    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[1,sl+1,c])/σ)))
end

for h=hh,sl=ssl-1, c=1:num_types, a=1:nwage
    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
end


for h=1:hh, c=1:num_types
    @addNLConstraint(m, q_est[h,c]==exp((EVA[h,1,c]-lease[h])/σ_s)/(exp((μ[work_time[h]])/σ_s)+exp((EVA[h,1,c]-lease[h])/σ_s)))
end;


for h=1:hh,sl=1:ssl-1, c=1:num_types
    @addNLConstraint(m, EVA[h,sl,c] == sum{pw[a]*EV[h,sl,c,a],a=1:nwage})
end


solve(m)