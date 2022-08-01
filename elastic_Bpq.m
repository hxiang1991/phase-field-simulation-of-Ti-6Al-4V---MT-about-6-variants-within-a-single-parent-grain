function [Bpq, Cijkl]= elastic_Bpq(e0,nx,ny,nz,kx,ky,kz,c11,c12,c44)

format long 

 %% == Assign Cijkl
 Cijkl(1:3, 1:3, 1:3, 1:3) = 0;
 
 Cijkl(1,1,1,1) = c11; Cijkl(2,2,2,2) = c11; Cijkl(3,3,3,3) = c11;
 Cijkl(1,1,2,2) = c12; Cijkl(2,2,3,3) = c12; Cijkl(3,3,1,1) = c12;
 Cijkl(2,2,1,1) = c12; Cijkl(3,3,2,2) = c12; Cijkl(1,1,3,3) = c12;
 
 Cijkl(1,2,1,2) = c44; Cijkl(2,3,2,3) = c44; Cijkl(3,1,3,1) = c44;
 Cijkl(2,1,2,1) = c44; Cijkl(3,2,3,2) = c44; Cijkl(1,3,1,3) = c44;
 Cijkl(1,2,2,1) = c44; Cijkl(3,2,2,3) = c44; Cijkl(1,3,3,1) = c44;
 Cijkl(2,1,1,2) = c44; Cijkl(2,3,3,2) = c44; Cijkl(3,1,1,3) = c44;

 %% == the corresponding stress tensor    
 s0 = zeros(3, 3, 6);
 for v= 1: 6
     for i = 1: 3
        for j = 1: 3
    
            s0(i,j,v) = Cijkl(i,j,1,1).* e0(1,1,v)+ Cijkl(i,j,2,2).* e0(2,2,v)+ Cijkl(i,j,3,3).* e0(3,3,v)+...
                     2*(Cijkl(i,j,1,2).* e0(1,2,v)+ Cijkl(i,j,1,3).* e0(1,3,v)+ Cijkl(i,j,2,3).* e0(2,3,v));
    
        end
     end
 end

 %% == iomega and omega
 iomega= zeros(3,3);
 Bpq= zeros(nx,ny,nz,6,6);

 for ix= 1: nx
     for iy= 1: ny
         for iz= 1: nx

             n1= kx(ix,iy,iz); n2= ky(ix,iy,iz); n3= kz(ix,iy,iz); n= [n1, n2, n3];

             for i = 1: 3
 	             for j = 1: 3
					               
	 	             iomega(i,j) = Cijkl(i,1,1,j).* n1.* n1+ Cijkl(i,1,2,j).* n1.* n2+ Cijkl(i,1,3,j).* n1.* n3+...
                                   Cijkl(i,2,1,j).* n2.* n1+ Cijkl(i,2,2,j).* n2.* n2+ Cijkl(i,2,3,j).* n2.* n3+...
                                   Cijkl(i,3,1,j).* n3.* n1+ Cijkl(i,3,2,j).* n3.* n2+ Cijkl(i,3,3,j).* n3.* n3;
                 end
             end

             omega = inv(iomega);

             for p = 1: 6
                 for q= 1: 6

                      bk0 = 0;
                      for i = 1: 3
                          for j = 1: 3
		                      for k = 1: 3
		 	                      for l = 1: 3
		 				 	                    
		 		                     bk0 = bk0+ Cijkl(i,j,k,l).* e0(i,j,p).* e0(k,l,q);
							                    
                                  end
                              end
                          end
                      end

                      Bpq(ix,iy,iz,p,q) = bk0 -n* s0(:,:,p)* omega* s0(:,:,q)* n';

                 end
             end

         end
     end
 end

 Bpq(isnan(Bpq))= 0; Bpq(isinf(Bpq))= 0;

end %  end function


