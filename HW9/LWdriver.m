function LWdriver(method)

mvals = [49 99 199];
N = length( mvals );

mesh_widths = zeros(1,N);
computed_errors = zeros(1,N);
for j = 1 : N
  if method=='lw' 
    [h,k,error] = advection_LW_pbc(mvals(j),N);
    mesh_widths(j) = h;
    computed_errors(j) = error;
  elseif method=='up'
    [h,k,error] = advection_up_pbc(mvals(j),N);
    mesh_widths(j) = h;
    computed_errors(j) = error;
end

error_table( mesh_widths, computed_errors );
error_loglog( mesh_widths, computed_errors );
end

