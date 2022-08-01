clear;clc
format long;

%% -- simulation paraters
nx= 128; ny= nx; nz= nx; nxyz= nx* ny* nz;      % simulation size
lengthx= nx* 2.5e-8;                            % length in x direction of simulation, unit: m
dx= lengthx/nx;                                 % mesh size, unit: m                                 
gamma= 50* 1e-3;                                % interfacial energy, unit: J/m2
Vmol= 1e-5;                                     % molar volume, unit: m3/mol
thick= 5* dx;                                   % thickness of interface, unit: m
ttime= 0;
Variants = 6;

dt0= 2e-03;  M0= 1.6e-07;                       % delta time and interfacial mobility
nstep= 2e04; nprint1= 1e02; nprint2= 1e03;      % for loop step and output step

diffGibbs= 3.1375e02;                           % driving force, unit: J/mol

kappa0= 3* gamma* thick/8^0.5;                  % gradient coefficient, unit: J/m
g_barrier= 3* Vmol* gamma/(4*2^0.5* thick);     % gibbs energy barrier, unit: J/mol

a0= 32* g_barrier;                              % gibbs free energy coefficient
b0= 3* a0+ 12* diffGibbs;
c0= 2* a0+ 12* diffGibbs;

ampnoise= 2e-3; noisestep= 680;                 % noise amptitude and nucleation stage

% -- dimensionless form
delta_Gm= diffGibbs; delta_M= 1e-18;
a= a0/delta_Gm; b= b0/delta_Gm; c= c0/delta_Gm;
kappa= kappa0* Vmol/(delta_Gm* dx^2);
dt= dt0* delta_Gm* delta_M/ dx^2;
M= M0* dx^2/(delta_M* Vmol);                    % interface mobility

% -- elastic strain tensor and elastic constants in parent orthonormal frame
U= zeros(3,3,6);                                % shape deformation
e0= zeros(3,3,6);                               % elastic strain

U(:,:,1) = [0.8925 0 0; 0 1.0474 -0.0456; 0 -0.0456 1.0474];  e0(:,:,1)= (U(:,:,1).'* U(:,:,1)- eye(3))/2;      
U(:,:,2) = [0.8925 0 0; 0 1.0474  0.0456; 0  0.0456 1.0474];  e0(:,:,2)= (U(:,:,2).'* U(:,:,2)- eye(3))/2; 
U(:,:,3) = [1.0474 0 -0.0456; 0 0.8925 0; -0.0456 0 1.0474];  e0(:,:,3)= (U(:,:,3).'* U(:,:,3)- eye(3))/2; 
U(:,:,4) = [1.0474 0  0.0456; 0 0.8925 0;  0.0456 0 1.0474];  e0(:,:,4)= (U(:,:,4).'* U(:,:,4)- eye(3))/2; 
U(:,:,5) = [1.0474 -0.0456 0; -0.0456 1.0474 0; 0 0 0.8925];  e0(:,:,5)= (U(:,:,5).'* U(:,:,5)- eye(3))/2; 
U(:,:,6) = [1.0474  0.0456 0;  0.0456 1.0474 0; 0 0 0.8925];  e0(:,:,6)= (U(:,:,6).'* U(:,:,6)- eye(3))/2; 

c11= 97.7e9/(delta_Gm/Vmol); c12= 82.7e9/(delta_Gm/Vmol); c44= 37.5e9/(delta_Gm/Vmol);

phiAplot1= zeros(nx, ny, nz); phiAplot2= zeros(nx, ny, nz); phiAplot3= zeros(nx, ny, nz); 
phiAplot4= zeros(nx, ny, nz); phiAplot5= zeros(nx, ny, nz); phiAplot6= zeros(nx, ny, nz); 
sumphiAplot= zeros(nx, ny, nz);

cc= hsv(Variants);
video= VideoWriter('microstructure evolution_twinning.avi');
video.FrameRate = 5; video.Quality= 100; open(video);

vflag=nstep/nprint1; VolF = zeros(vflag+ 1, Variants+ 1);
VolF(2:vflag+ 1,1)= (1: vflag)* nprint1;

[sym_x, sym_y, sym_z] = ndgrid(1: nx, 1: ny, 1: nz);
tmp_x = reshape(sym_x, nx* ny* nz, 1);
tmp_y = reshape(sym_y, nx* ny* nz, 1);
tmp_z = reshape(sym_z, nx* ny* nz, 1);
sym_cor_mat = [tmp_x tmp_y tmp_z];

phiA= zeros(nx, ny, nz, Variants);

tmpkx= 2*pi*[0: nx/2 -nx/2+1: -1]/nx; tmpky= tmpkx; tmpkz= tmpkx;
[kx,ky,kz]= ndgrid(tmpkx,tmpky,tmpkz);k2= kx.^2+ ky.^2+ kz.^2;
kx= kx./k2.^0.5; ky= ky./k2.^0.5; kz= kz./k2.^0.5;
kx(isnan(kx))= 0; ky(isnan(ky))= 0; kz(isnan(kz))= 0; 

[Bpq, Cijkl]= elastic_Bpq(e0,nx,ny,nz,kx,ky,kz,c11,c12,c44);

denom= 1+ dt* M* kappa* k2;

for istep= 1: nstep

	ttime= ttime+ dt;

    phiA2= phiA.^2; sumphiA2= sum(phiA2, 4);
    
    phiAk= fft(fft(fft(phiA, [], 1), [], 2), [], 3);

    dfdphiA= a* phiA- b* phiA2+ c* phiA.* sumphiA2;
    dfdphiAk= fft(fft(fft(dfdphiA, [], 1), [], 2), [], 3);
    
    deldphiAk(:,:,:,1)= Bpq(:,:,:,1,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,1,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,1,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,1,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,1,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,1,6).* phiAk(:,:,:,6); 
    deldphiAk(:,:,:,2)= Bpq(:,:,:,2,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,2,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,2,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,2,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,2,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,2,6).* phiAk(:,:,:,6); 
    deldphiAk(:,:,:,3)= Bpq(:,:,:,3,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,3,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,3,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,3,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,3,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,3,6).* phiAk(:,:,:,6); 
    deldphiAk(:,:,:,4)= Bpq(:,:,:,4,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,4,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,4,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,4,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,4,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,4,6).* phiAk(:,:,:,6); 
    deldphiAk(:,:,:,5)= Bpq(:,:,:,5,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,5,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,5,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,5,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,5,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,5,6).* phiAk(:,:,:,6); 
    deldphiAk(:,:,:,6)= Bpq(:,:,:,6,1).* phiAk(:,:,:,1)+ Bpq(:,:,:,6,2).* phiAk(:,:,:,2)+ Bpq(:,:,:,6,3).* phiAk(:,:,:,3)+...
                        Bpq(:,:,:,6,4).* phiAk(:,:,:,4)+ Bpq(:,:,:,6,5).* phiAk(:,:,:,5)+ Bpq(:,:,:,6,6).* phiAk(:,:,:,6); 

    if istep< noisestep
        
        r1= rand(nx,ny,nz,Variants); r2= rand(nx,ny,nz,Variants);
        noise= -2* log(r1).* sin(2*pi*r2).*cos(2*pi*r2)*ampnoise;
        noisek= fft(fft(fft(noise, [], 1), [], 2), [], 3);
		
        phiAk= (phiAk- dt* M.* (dfdphiAk+ deldphiAk))./denom+ noisek;

    else

        phiAk= (phiAk- dt* M.* (dfdphiAk+ deldphiAk))./denom;

    end

	phiA= real(ifft(ifft(ifft(phiAk, [], 1), [], 2), [], 3));

    inrange = (phiA>1); phiA(inrange)= 1;
    inrange = (phiA<0); phiA(inrange)= 0;
    
    phiA2= phiA.^2;
    tmpphiA1= phiA2(:,:,:,1); tmpphiA2= phiA2(:,:,:,2); tmpphiA3= phiA2(:,:,:,3); 
    tmpphiA4= phiA2(:,:,:,4); tmpphiA5= phiA2(:,:,:,5); tmpphiA6= phiA2(:,:,:,6); 

    inrange= (tmpphiA1>0.5); phiAplot1(inrange)= 1; inrange= (tmpphiA2>0.5); phiAplot2(inrange)= 1; inrange= (tmpphiA3>0.5); phiAplot3(inrange)= 1;
    inrange= (tmpphiA4>0.5); phiAplot4(inrange)= 1; inrange= (tmpphiA5>0.5); phiAplot5(inrange)= 1; inrange= (tmpphiA6>0.5); phiAplot6(inrange)= 1;
    sumphiAplot= phiAplot1+ phiAplot2+ phiAplot3+ phiAplot4+ phiAplot5+ phiAplot6;
    
    if (mod(istep, nprint1)== 0)
         
        f= figure('visible','off');
        fv1= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,1), 0.5); patch(fv1, 'FaceColor', cc(1, 1: 3), 'EdgeColor', 'none'); hold on  
        fv2= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,2), 0.5); patch(fv2, 'FaceColor', cc(2, 1: 3), 'EdgeColor', 'none'); hold on  
        fv3= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,3), 0.5); patch(fv3, 'FaceColor', cc(3, 1: 3), 'EdgeColor', 'none'); hold on  
        fv4= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,4), 0.5); patch(fv4, 'FaceColor', cc(4, 1: 3), 'EdgeColor', 'none'); hold on 
        fv5= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,5), 0.5); patch(fv5, 'FaceColor', cc(5, 1: 3), 'EdgeColor', 'none'); hold on  
        fv6= isosurface(sym_x, sym_y, sym_z, phiA2(:,:,:,6), 0.5); patch(fv6, 'FaceColor', cc(6, 1: 3), 'EdgeColor', 'none'); hold on 
        axis equal; axis([1, nx, 1, ny, 1, nz]); grid on; view(3);
        box on; ax= gca; ax.BoxStyle= 'full'; xlabel('x'); ylabel('y');zlabel('z');
        ax.XColor= [0.9, 0.9, 0.9]; ax.YColor= [0.9, 0.9, 0.9]; ax.ZColor = [0.9, 0.9, 0.9];
        set(gca,'color',[0.3 0.3 0.3]);set(gcf,'color',[0.3 0.3 0.3]);
        camlight; lighting phong;title(['time: ', num2str(ttime,'%6.3f'), 's'], 'Color', 'w');
        
        legend({'V1', 'V2', 'V3', 'V4', 'V5', 'V6'},'NumColumns',1, 'TextColor','w'); legend boxoff; 
		
		frame= getframe(gcf); writeVideo(video, frame); 
		
		VolF(istep/nprint1+ 1,2)= sum(phiAplot1, 'all')* 100/nxyz; VolF(istep/nprint1+ 1,3)= sum(phiAplot2, 'all')* 100/nxyz; 
		VolF(istep/nprint1+ 1,4)= sum(phiAplot3, 'all')* 100/nxyz; VolF(istep/nprint1+ 1,5)= sum(phiAplot4, 'all')* 100/nxyz; 
		VolF(istep/nprint1+ 1,6)= sum(phiAplot5, 'all')* 100/nxyz; VolF(istep/nprint1+ 1,7)= sum(phiAplot6, 'all')* 100/nxyz;
		
		if (mod(istep, nprint2)== 0)
		
			filename= sprintf('time_%6.3f.fig', ttime);  savefig(f, filename); 
     
			phiA2_1= phiA2(:,:,:,1); phiA2_2= phiA2(:,:,:,2); phiA2_3= phiA2(:,:,:,3);
			phiA2_4= phiA2(:,:,:,4); phiA2_5= phiA2(:,:,:,5); phiA2_6= phiA2(:,:,:,6); 
			phiAplot= 1* phiAplot1+ 2* phiAplot2+ 3* phiAplot3+ 4* phiAplot4+ 5* phiAplot5+ 6* phiAplot6;
		
			write_vtk_grid_values(nx,ny,nz,sym_cor_mat,dx,ttime,phiAplot,phiA2_1,phiA2_2,phiA2_3,phiA2_4,phiA2_5,phiA2_6);
	
		end
	
		f;clf;
    
    end
	
	

    sumphiA2= sum(phiA2, 4);max_phiA2= max(sumphiA2,[],'all'); 
    maxphiAplot= max(sumphiAplot,[],'all');
              
    phiAplot1= zeros(nx, ny, nz); phiAplot2= zeros(nx, ny, nz); phiAplot3= zeros(nx, ny, nz);
    phiAplot4= zeros(nx, ny, nz); phiAplot5= zeros(nx, ny, nz); phiAplot6= zeros(nx, ny, nz); 

    if (maxphiAplot> 1 |max_phiA2< 1e-4| isinf(deldphiAk))				
		break				
    end
end

filename= 'VolumnFraction.mat';
save(filename, 'VolF'); close(video);