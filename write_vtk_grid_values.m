function [  ]= write_vtk_grid_values(nx,ny,nz,sym_cor_mat,dx,ttime,data1,data2,data3,data4,data5,data6,data7)         

format short;

%% == open output file

fname = sprintf('time_%6.3f.vtk', ttime);
out = fopen(fname, 'w');

npoint = nx* ny* nz;

%% == reshape the Nx-Ny-Nz matrix into Nx* Ny* Nz-1 array for convenient output   

data1 = reshape(data1, nx* ny* nz, 1);
data2 = reshape(data2, nx* ny* nz, 1); data3 = reshape(data3, nx* ny* nz, 1);
data4 = reshape(data4, nx* ny* nz, 1); data5 = reshape(data5, nx* ny* nz, 1);
data6 = reshape(data6, nx* ny* nz, 1); data7 = reshape(data7, nx* ny* nz, 1);

%% == start writing ASCII VTK file:

% == header of VTK file

fprintf(out, '# vtk DataFile Version 2.0\n');
fprintf(out, 'time_10.vtk\n');
fprintf(out, 'ASCII\n');
fprintf(out, 'DATASET STRUCTURED_GRID\n');

% == coords of grid points:

fprintf(out, 'DIMENSIONS %5d  %5d  %5d\n', nx, ny, nz);
fprintf(out, 'POINTS %7d   float\n', npoint);

sym_cor_mat = sym_cor_mat* dx;
Array_sym_cor_mat = reshape(sym_cor_mat.', nx* ny* nz* 3, 1);

formatSpec = '%14.6e  %14.6e  %14.6e\n';
fprintf(out, formatSpec, Array_sym_cor_mat);

% == write grid point values:

fprintf(out,'POINT_DATA %5d\n',npoint);

fprintf(out,'SCALARS variants_index  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data1);                % output grid value of variants index
            
fprintf(out,'SCALARS variant_1  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data2);                 % output grid value of variant 1

fprintf(out,'SCALARS variant_2  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data3);                 % output grid value of variant 2
            
fprintf(out,'SCALARS variant_3  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data4);                 % output grid value of variant 3
            
fprintf(out,'SCALARS variant_4  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data5);                 % output grid value of variant 4
            
fprintf(out,'SCALARS variant_5  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data6);                 % output grid value of variant 5

fprintf(out,'SCALARS variant_6  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data7);                 % output grid value of variant 6

fclose(out);

end %endfunction

      
