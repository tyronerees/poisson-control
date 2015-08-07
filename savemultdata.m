% multdata = genmultdata(10,'dirichlet');
% save multdataQ1_dir  multdata
% multdata = genmultdata(10,'neumann');
% save multdataQ1_neu  multdata
% multdata = genmultdata(10,'mixed');
% save multdataQ1_mix multdata
multdata = genmultdata(10,'neumannbc');
save multdataQ1_neubc  multdata
% pow = 6;
% N = 2^(pow);
% cd PoissonControl
% [x y z cub bdy] = mesh3d_bdy(N);
% cd ..
% multdata = genmultdata3d(pow,unique(bdy'));
% save multdataQ1_dir_3d  multdata
