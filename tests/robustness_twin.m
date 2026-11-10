%ROBUSTNESS_TWIN  Robustness of the CG ranking to the noise level (Table 3).
%   Runs FR, PR, PB and DY on synthetic data with known D and S at three noise
%   levels and prints the decrease of J, the errors of D and S and the number
%   of restarts (same grid and fields as tests/test_twin_synthetic.m).
%   Uses the functions in ../core (added to the path below); base MATLAB only.
%   Part of ychannel-adjoint-inverse; see README.md for the article reference.
clear; clc; here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'..','core')); addpath(here);
Nx=81; Ny=270; x=linspace(-0.5,0.5,Nx); y=(0:Ny-1).'*0.03; [X,Y]=meshgrid(x,y); v=0.15;
Db=[1.76e-3,2.03e-3,0.81e-3]; Gd=exp(-X.^2/(2*0.03^2)); Gs=exp(-X.^2/(2*0.025^2)); a=0.05*exp(-Y/4);
Dt={Db(1)*(1+1.5*Gd),Db(2)*(1+1.5*Gd),Db(3)*(1+1.5*Gd)}; St={-a.*Gs,-2*a.*Gs,a.*Gs};
Cin={1.0*0.5*erfc(X/0.01),2.0*0.5*erfc(-X/0.01),zeros(Ny,Nx)}; W={ones(Ny,Nx),ones(Ny,Nx),ones(Ny,Nx)};
opt=struct('alphaD',1e-6,'alphaS',1e-6,'Dref',1e-3,'Sref',1e-2,'Dmin',1e-5,'Dmax',1e-1);
p0=make_problem(x,y,v,Cin,W,Cin,opt); Ct=cell(1,3); for s=1:3, Ct{s}=forward_solve(Dt{s},St{s},s,p0); end
D0={Db(1)*ones(Ny,Nx),Db(2)*ones(Ny,Nx),Db(3)*ones(Ny,Nx)}; S0=cell(1,3);
for s=1:3, S0{s}=init_source_balance(Ct{s},x,y,v,D0{s},0.03); end
band=abs(X)<0.08 & Y>0.5; relerr=@(A,B) norm(A(band)-B(band))/norm(B(band));
methods={'FR','PR','PB','DY'}; noises=[0.01 0.03 0.087];
fprintf('%6s %4s %8s %11s %11s %7s %7s %8s\n','noise','CG','drop%','J_end','Jdata','errD','errS','restarts');
for nz=noises
  Cobs=cell(1,3); for s=1:3, Cobs{s}=Ct{s}+nz*det_randn([Ny Nx],10+s); end
  prob=make_problem(x,y,v,Cobs,W,Ct,opt);
  for m=1:4
    r=cg_optimize(prob,D0,S0,methods{m},struct('maxIter',100,'verbose',false));
    eD=mean([relerr(r.D1,Dt{1}),relerr(r.D2,Dt{2}),relerr(r.D3,Dt{3})]);
    eS=mean([relerr(r.S1,St{1}),relerr(r.S2,St{2}),relerr(r.S3,St{3})]);
    fprintf('%6.3f %4s %8.2f %11.4e %11.4e %7.3f %7.3f %8d\n',nz,methods{m},100*(1-r.final_cost/r.initial_cost),r.final_cost,r.J_data_history(end),eD,eS,r.cg_restart_count);
  end
end
