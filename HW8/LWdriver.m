function LWdriver(method)

mvals = [199 399 799 1999 3999];
N = length( mvals );

mesh_widths = zeros(1,N);
computed_errors = zeros(1,N);
for j = 1 : N
  if method=='lfa' 
    [h,k,error] = advection_lf_pbc(mvals(j));
    mesh_widths(j) = h;
    computed_errors(j) = error;
  elseif method=='lfb'
    [h,k,error] = advection_lf_pbc_b(mvals(j));
    mesh_widths(j) = h;
    computed_errors(j) = error;
  end
end

if method=='lfc'
    m = mvals(1);
    [h,k,error,xvec,tvec,domega] = advection_lf_pbc_c(m);
    mesh_widths(j) = h;
    computed_errors(j) = error;
    v = polyfit(tvec,xvec,1);
    disp(domega)
    disp(v)
end
error_table( mesh_widths, computed_errors );
error_loglog( mesh_widths, computed_errors );
end