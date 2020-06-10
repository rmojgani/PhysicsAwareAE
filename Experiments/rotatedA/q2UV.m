function [ U_d , V_d , d] = q2UV( q, Nx, Nt, kk, Flag)
% Extracts and reshape basis functions vectorized in vector q based on its
% size

U_d = reshape(q(            1:Flag(1)                   ,1),Nx,kk);
V_d = reshape(q(           Flag(1)+1:sum(Flag(1:2))     ,1),kk,Nt);

d = U_d*V_d;
end
%% Attemp to generalizing the function ... yet to be ... 
% Flag=[0,Flag] % Add a zero for sake of the loop
% for icount = 1:1:numel(Flag)
%      U.(['v_',num2str(1)]) = reshape(q(           sum(Flag(1:icount))+1:Flag(icount+1)                   ,1),Nx*Ny,k);
%      V.(['v_',num2str(1)]) = reshape(q(            1:Flag(1)                   ,1),k,Nx*Ny);

% end
% U_w = reshape(q(            1:Flag(1)                   ,1),Nx*Ny,k);
% V_w = reshape(q(           Flag(1)+1:sum(Flag(1:2))     ,1),k,Nt);
% 
% U_x = reshape(q(    sum(Flag(1:2))+1:sum(Flag(1:3))     ,1),Nx*Ny,k);
% V_xy = reshape(q(    sum(Flag(1:3))+1:sum(Flag(1:4))     ,1),k,Nt);
% 
% U_y = reshape(q(    sum(Flag(1:4))+1:sum(Flag(1:5))     ,1),Nx*Ny,k);