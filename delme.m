
sz_global_idx_map = 2;
f0u=zeros(18,18);
%f0u=zeros(24,24)
 n=1;
%  for j=1:sz_global_idx_map
%      %tmp=reshape(1:81,9,9);
%      tmp=reshape(1:64,8,8);
%      sz_tmp = size(tmp,1);
%      m=0;
%      for k = 1:sz_global_idx_map
%          j = n;
%          f0u(j:sz_tmp+n-1, m+1:m+sz_tmp ) = tmp;
%          m = m+sz_tmp;
%      end
%      n = n+sz_tmp;
%  end
 
        n=1;
         s=1;
         for j=1:sz_global_idx_map
             m=0;
             for k = 1:sz_global_idx_map
                 %tmp = reshape(1:64,8,8);
                 tmp=reshape(1:81,9,9);
                 sz_tmp = size(tmp,1);
                 j = n;
                 f0u(j:sz_tmp+n-1, m+1:m+sz_tmp ) = tmp;
                 m = m+sz_tmp;
                 s =s+1;
             end
             n = n+sz_tmp;
         end