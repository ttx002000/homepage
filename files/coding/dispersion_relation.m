clear all;
hold off;
N1=40;
k_path=zeros(3,N1);
m=23;
n=22;%make sure that m>n
for ja=1:N1
k_path(:,ja)=(m*[1;0;0]+n*[1/2;3^(1/2)/2;0])*ja/(N1-1)*4*pi/(3*(m^2+n^2+m*n));
end
%for ja=1:N1
%k_path(:,ja)=((m-n)*[1;0;0]+(m+2*n)*[1/2;3^(1/2)/2;0])*ja/(N1-1)*(1/((m-n)^2+(m+2*n)^2+(m-n)*(m+2*n)))^(1/2)*2*pi/(3*(m^2+n^2+m*n))^(1/2);
%end
k_length=zeros(1,N1);

for ja=1:N1
k_length(ja)=dot(k_path(:,ja),k_path(:,ja))^(1/2);
end

theta=acos(1/2*(m^2+n^2+4*m*n)/(m^2+n^2+m*n));
V_pi=-2.7;
V_sg=0.48;
a_1=[1;0;0];
a_2=[1/2;3^(1/2)/2;0];
d=0.335/0.246;

band=cell(1,N1);

num_atom=0;
col_pos=zeros(3,16*(m+n)^2);

for ma=-(m+n):(m+n)
    for na=0:2*(m+n)
  
   
   if ((ma*(m+n)+na*n)<(m^2+n^2+m*n)) && ((ma*(m+n)+na*n)>=0) ...
       && ((na*m-ma*n)<(m^2+n^2+m*n)) && ((na*m-ma*n)>=0)
   
   num_atom=num_atom+1;
   col_pos(:,num_atom)=ma*[1;0;0]+na*[1/2;3^(1/2)/2;0];
   end
   
    end
end

col_pos(:,num_atom+1:num_atom*2)=col_pos(:,1:num_atom)+[1/2;1/(2*3^(1/2));0];
num_atom=num_atom*2;



num_count=0;
for ma=-(m+n):(m+n)
    for na=0:2*(m+n)
   if ((ma*(m+n)+na*m)/(m^2+n^2+m*n)<1) && ((ma*(m+n)+na*m)/(m^2+n^2+m*n)>=0 )...
       && ((na*n-ma*m)/(m^2+n^2+m*n)<1) && ((na*n-ma*m)/(m^2+n^2+m*n)>=0)
   num_atom=num_atom+1;
   col_pos(:,num_atom)=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1]*(ma*[1;0;0]+na*[1/2;3^(1/2)/2;0]+[0;0;d]);
  num_count=num_count+1;
   end
   
    end
end

col_pos(:,num_atom+1:num_atom+num_count)=col_pos(:,num_atom-num_count+1:num_atom)+ ...
[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1]*[1/2;1/(2*3^(1/2));0];

num_atom=num_atom+num_count;

col_pos(:,num_atom+1:16*(m+n)^2)=[];


count=0;
row_pos=zeros(3,9*num_atom);

for ja=-1:1
    for jb=-1:1
        
   row_pos(:,1+num_atom*count:num_atom+num_atom*count)=col_pos(:,1:num_atom)+...
       ja*(m*[1;0;0]+n*[1/2;3^(1/2)/2;0])+jb*(-n*[1;0;0]+(n+m)*[1/2;3^(1/2)/2;0]);
    
    count=count+1;
    end
end



parfor ja=1:N1
    
  
   
    
       
    
    H_kspace=zeros(num_atom,num_atom);
    
       for jb=1:num_atom
           for jc=1:num_atom
               for jd=0:8
              distance=dot(row_pos(:,jc+jd*num_atom)-col_pos(:,jb),row_pos(:,jc+jd*num_atom)-col_pos(:,jb))^(1/2);
            
           if (distance<=4/3^(1/2)) &&  (distance~=0)%here we must be careful with onsite energy
              V_pi_sp=V_pi*exp(-(distance-1/3^(1/2))/0.184);
              V_sg_sp=V_sg*exp(-(distance-d)/0.184);
            H_kspace(jb,jc)=H_kspace(jb,jc)+(V_pi_sp*(1-dot(row_pos(:,jc+jd*num_atom)-col_pos(:,jb),[0;0;1])^2/distance^2)...
            +V_sg_sp*dot(row_pos(:,jc+jd*num_atom)-col_pos(:,jb),[0;0;1])^2/distance^2)...
            *exp(1i*dot(k_path(:,ja),row_pos(:,jc+jd*num_atom)-col_pos(:,jb)));
           end
               end
               end
       end
       
       
       A=eig(H_kspace);
     A=real(A);
     band{ja}=sort(A);
     
    
 
     
    
end
eigenvalue=zeros(N1,num_atom);

for ja=1:N1
eigenvalue(ja,:)=band{ja};
end

for ja=1:num_atom
plot(k_length,real(eigenvalue(:,ja)));
hold on;
end
ylim([-1 1])



