function D=FD(N,m,k,h)

Ax=zeros(m,m);
Au=zeros(m,m);
Av=zeros(m,m);
D=zeros(N,N);
b=zeros(m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m==2 % 2nd order approx
   b(k)=1;
   for i=1:2
      for j=1:2        
         if j>1            
            Ax(i,j)=(j-1).^i;            
         else            
            Ax(i,j)=(j-2).^i;                              
         end  
         Au(i,j)=j.^i;
         Av(i,j)=(j-3).^i;
      end            
   end
   % interior
   x=inv(Ax)*b; 
   sum1=sum(x); sum2=sum( Ax(k,:)*x ); x=x./sum2;
   for i=2:N-1
      D(i,i)=-sum1./sum2; 
      D(i,i-1)=x(1);
      D(i,i+1)=x(2);
   end
   
   u=inv(Au)*b; 
   sum1=sum(u); sum2=sum( Au(k,:)*u ); u=u./sum2;
   D(1,1)=-sum1./sum2; D(1,2)=u(1); D(1,3)=u(2);
   
   v=inv(Av)*b;    
   sum1=sum(v); sum2=sum( Av(k,:)*v ); v=v./sum2;
   D(N,N-2)=v(1); D(N,N-1)=v(2); D(N,N)=-sum1./sum2;   
   
end      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if m==4 % 4th order approx
   b(k)=1;
   for i=1:4
      for j=1:4        
         if j>2            
            Ax(i,j)=(j-2).^i;            
         else            
            Ax(i,j)=(j-3).^i;                              
         end         
         Au1(i,j)=j.^i;
         if j>1
            Au2(i,j)=(j-1).^i;
         else
            Au2(i,j)=(j-2).^i;
         end         
         Av1(i,j)=(j-5).^i;
         if j>3
            Av2(i,j)=(j-3).^i;
         else
            Av2(i,j)=(j-4).^i;
         end                  
      end         
   end
   % interior
   x=inv(Ax)*b; 
   sum1=sum(x); sum2=sum( Ax(k,:)*x ); x=x./sum2;
   for i=3:N-2
      D(i,i)=-sum1./sum2; 
      D(i,i-2)=x(1);
      D(i,i-1)=x(2);
      D(i,i+1)=x(3);
      D(i,i+2)=x(4);   
   end
   % boundary
   u1=inv(Au1)*b; 
   sum1=sum(u1); sum2=sum( Au1(k,:)*u1 ); u1=u1./sum2;  
   D(1,1)=-sum1./sum2; 
   D(1,2)=u1(1); D(1,3)=u1(2);
   D(1,4)=u1(3); D(1,5)=u1(4);
   u2=inv(Au2)*b; 
   sum1=sum(u2); sum2=sum( Au2(k,:)*u2 ); u2=u2./sum2;  
   D(2,1)=u2(1); 
   D(2,2)=-sum1./sum2; D(2,3)=u2(2);
   D(2,4)=u2(3); D(2,5)=u2(4);
   
   v1=inv(Av1)*b;    
   sum1=sum(v1); sum2=sum( Av1(k,:)*v1 ); v1=v1./sum2;
   D(N,N-4)=v1(1); D(N,N-3)=v1(2); 
   D(N,N-2)=v1(3); D(N,N-1)=v1(4);    
   D(N,N)=-sum1./sum2;   
   v2=inv(Av2)*b;    
   sum1=sum(v2); sum2=sum( Av2(k,:)*v2 ); v2=v2./sum2;
   D(N-1,N-4)=v2(1); D(N-1,N-3)=v2(2); 
   D(N-1,N-2)=v2(3); D(N-1,N-1)=-sum1./sum2;    
   D(N-1,N)=v2(4);   
end      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m==6 % 6th order approx
   b(k)=1;
   for i=1:6
      for j=1:6        
         if j>3            
            Ax(i,j)=(j-3).^i;            
         else            
            Ax(i,j)=(j-4).^i;                              
         end         
         Au1(i,j)=j.^i;
         if j>1
            Au2(i,j)=(j-1).^i;
         else
            Au2(i,j)=(j-2).^i;
         end 
         if j>2
            Au3(i,j)=(j-2).^i;
         else
            Au3(i,j)=(j-3).^i;
         end          
         Av1(i,j)=(j-7).^i;
         if j>5
            Av2(i,j)=(j-5).^i;
         else
            Av2(i,j)=(j-6).^i;
         end                  
         if j>4
            Av3(i,j)=(j-4).^i;
         else
            Av3(i,j)=(j-5).^i;
         end                           
      end         
   end
   % interior
   x=inv(Ax)*b; 
   sum1=sum(x); sum2=sum( Ax(k,:)*x ); x=x./sum2;
   for i=4:N-3
      D(i,i)=-sum1./sum2; 
      D(i,i-3)=x(1);      
      D(i,i-2)=x(2);
      D(i,i-1)=x(3);
      D(i,i+1)=x(4);
      D(i,i+2)=x(5);   
      D(i,i+3)=x(6);      
   end
   % boundary
   u1=inv(Au1)*b; 
   sum1=sum(u1); sum2=sum( Au1(k,:)*u1 ); u1=u1./sum2;  
   D(1,1)=-sum1./sum2; 
   D(1,2)=u1(1); D(1,3)=u1(2);
   D(1,4)=u1(3); D(1,5)=u1(4);
   D(1,6)=u1(5); D(1,7)=u1(6);   
   u2=inv(Au2)*b; 
   sum1=sum(u2); sum2=sum( Au2(k,:)*u2 ); u2=u2./sum2;  
   D(2,1)=u2(1); 
   D(2,2)=-sum1./sum2; D(2,3)=u2(2);
   D(2,4)=u2(3); D(2,5)=u2(4);
   D(2,6)=u2(5); D(2,7)=u2(6);   
   u3=inv(Au3)*b; 
   sum1=sum(u3); sum2=sum( Au3(k,:)*u3 ); u3=u3./sum2;  
   D(3,1)=u3(1); 
   D(3,2)=u3(2); D(3,3)=-sum1./sum2;
   D(3,4)=u3(3); D(3,5)=u3(4);
   D(3,6)=u3(5); D(3,7)=u3(6);   
   v1=inv(Av1)*b;    
   sum1=sum(v1); sum2=sum( Av1(k,:)*v1 ); v1=v1./sum2;
   D(N,N-6)=v1(1); D(N,N-5)=v1(2);    
   D(N,N-4)=v1(3); D(N,N-3)=v1(4); 
   D(N,N-2)=v1(5); D(N,N-1)=v1(6);    
   D(N,N)=-sum1./sum2;
   v2=inv(Av2)*b;    
   sum1=sum(v2); sum2=sum( Av2(k,:)*v2 ); v2=v2./sum2;
   D(N-1,N-6)=v2(1); D(N-1,N-5)=v2(2);    
   D(N-1,N-4)=v2(3); D(N-1,N-3)=v2(4); 
   D(N-1,N-2)=v2(5); D(N-1,N-1)=-sum1./sum2 ;    
   D(N-1,N)=v2(6);   
   v3=inv(Av3)*b;    
   sum1=sum(v3); sum2=sum( Av3(k,:)*v3 ); v3=v3./sum2;
   D(N-2,N-6)=v3(1); D(N-2,N-5)=v3(2);    
   D(N-2,N-4)=v3(3); D(N-2,N-3)=v3(4); 
   D(N-2,N-2)=-sum1./sum2; D(N-2,N-1)=v3(5);    
   D(N-2,N)=v3(6);   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const=(factorial(k)./ (h.^k) );
D=const.*D;      
   




