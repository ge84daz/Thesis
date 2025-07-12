function [] = updateWaitbar(Nx_hs,wait,x,Nf)

global numCompleted

numCompleted = numCompleted + size(x,1);
waitbar(numCompleted/Nx_hs/Nf,wait,sprintf('Progress: %1.1f %%',numCompleted/Nx_hs/Nf*100));

end

