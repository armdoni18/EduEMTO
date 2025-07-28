function [Bx, By, B, MagEn, TotalMagEn] = Func6_MagFlux_BxBy(fem)

    Bx      = zeros(fem.ne,1);
    By      = zeros(fem.ne,1);
    B       = zeros(fem.ne,1);
    MagEn   = zeros(fem.ne,1);

    for e=1:fem.ne
        i = fem.IX(e,1);  j = fem.IX(e,2);  k = fem.IX(e,3);

        xi = fem.X(i,1);  yi = fem.X(i,2);
        xj = fem.X(j,1);  yj = fem.X(j,2);
        xk = fem.X(k,1);  yk = fem.X(k,2);

        bi = yj - yk;  ci = xk - xj;
        bj = yk - yi;  cj = xi - xk;
        bk = yi - yj;  ck = xj - xi;

        % Compute magnetic flux density
        Bx(e)   = ( 1 / (2 * fem.Ve(e))) * (ci * fem.A(i) + cj * fem.A(j) + ck * fem.A(k));
        By(e)   = (-1 / (2 * fem.Ve(e))) * (bi * fem.A(i) + bj * fem.A(j) + bk * fem.A(k));
        B(e)    = sqrt(Bx(e)^2 + By(e)^2);

        % Magnetic energy computation
        mu_e    = fem.IX(e,5);  
        MagEn(e)= B(e)^2 * fem.Ve(e) * mu_e/2;
    end

        % Store to fem struct
        fem.Bx  = Bx;
        fem.By  = By;
        fem.B   = B;

        fem.MagEn       = MagEn;
        TotalMagEn      = sum(MagEn);
        fem.TotalMagEn  = TotalMagEn; 

    fprintf('Magnetic flux density computation Done. ✅\n');
end
