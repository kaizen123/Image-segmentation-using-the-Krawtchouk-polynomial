% CREASEG_SUVADIP: Version 1 of Smooth Basis Chan Vese 
% Modified by Kambiz Rahbar, 2016.
%------------------------------------------------------------------------
% Legendre Basis Chan Vese: Faster implementation with vectorized basis
% Inputs: Featured_Img 2D image
%         mask        Initialization (1 = foreground, 0 = bg)
%         max_its     Number of iterations to run segmentation for
%         length_term (optional) Weight of smoothing term; higer = smoother.  default = 0.2
%         poly_degree Degree of the Legendre basis polynomial  
%         display     (optional) displays intermediate outputs; default = true
% Outputs: seg        Final segmentation mask (1=fg, 0=bg)
%
% Description: The internal and external intensities are modeled as linear combination of legendre basis function, making it more adaptive to
% intensity inhomogeneity. Furthermore, the basis coefficients can be obtained by a close form solution, which makes it computationally
% feasible. This is the computationally faster version, where the basis vectors are vectorized.
% -------------------------------------------------------------------------
% IF YOU ARE USING THIS CODE, PLEASE CITE THE FOLLOWING PAPER:
% Mukherjee, S.; Acton, S., "Region based segmentation in presence of intensity 
% inhomogeneity using Legendre polynomials," Signal Processing Letters, IEEE , vol.PP, no.99, pp.1,1 doi: 10.1109/LSP.2014.2346538

%% Krawtchouk
function [seg, evolved_phi, its] = creaseg_Krawtchouk(Img, init_mask, prm)
    prm.init_phi = computeSDF(init_mask); % phi: weighted mask; inside object>0, outside<0, border=0
    prm.basis_vec = Gen_2D_Vec_Basis_Fucn(prm.Featured_Img, prm.poly_degree, prm.display); % each col contains a basis vector (Legendre, ...) of length = vec(Featured_Img)
    prm.Img = Img;
    prm.display_intrvl = 30;
    prm.lambda1 = 1;
    prm.lambda2 = 1;
    prm.convg_count = 20;

    [evolved_phi, ~, ~, its] = ChanVeseCustFunc_vectorized(prm);

    seg = evolved_phi >= 0;
    evolved_phi = -evolved_phi;     % Return the negative following Creaseg's convention
    
    % by kambiz rahbar
    save_resImage('krawtchouk');
  
end

%% ChanVeseCustFunc_vectorized
function [phi, masked_Fg_est_obj, masked_Bg_est_obj, its] = ChanVeseCustFunc_vectorized(prm)
    phi_update_rate = 1;
    
    Heaviside = @(y,e) (0.5 * (1+ (2/pi) * atan(y/e) ) );
    Dirac_global = @(x,e) ((e/pi)./(e^2.+ x.^2)); % Dirac_global is defrencial of Heaviside function

    phi = prm.init_phi;

    stop = 0;
    count = 0; % count how many succesive times we have attained convergence, reduce local minima
    II = eye(size(prm.basis_vec, 2), size(prm.basis_vec, 2));
    its = 0;
    while (its < prm.max_iter && stop == 0)

        h_phi = Heaviside(phi,2); % the MASK, 
        Fg_mask = h_phi;   % foreground object(s) mask
        Bg_mask = 1-h_phi; % background object(s) mask

        Fg_Featured_Img = prm.Featured_Img .* Fg_mask; % image mased by foreground object(s) mask
        Bg_Featured_Img = prm.Featured_Img .* Bg_mask; % image mased by background object(s) mask

        Fg_basis_vec = prm.basis_vec .* repmat(Fg_mask(:), 1, size(prm.basis_vec,2)); % basis vector masked by foreground object(s) mask
        Bg_basis_vec = prm.basis_vec .* repmat(Bg_mask(:), 1, size(prm.basis_vec,2)); % basis vector masked by background object(s) mask
        
        % the foreground and background needed to be estimated by foreground and background basis vectors.
        % the parameters are:
        Fg_est_prm = pinv(prm.basis_vec' * Fg_basis_vec + prm.l2_reg_term * II) * prm.basis_vec' * Fg_Featured_Img(:);
        Bg_est_prm = pinv(prm.basis_vec' * Bg_basis_vec + prm.l2_reg_term * II) * prm.basis_vec' * Bg_Featured_Img(:);

        Fg_est_obj = prm.basis_vec * Fg_est_prm; % Foreground object estimated by the use of basis_vec (A'*p(x) in eq.3)
        Bg_est_obj = prm.basis_vec * Bg_est_prm; % Background object estimated by the use of basis_vec

        Fg_est_obj = reshape(Fg_est_obj, size(prm.Featured_Img));
        Bg_est_obj = reshape(Bg_est_obj, size(prm.Featured_Img));

        masked_Fg_est_obj = Fg_est_obj .* Fg_mask; % Foreground estimated object masked by foreground mask (Func return val.)
        masked_Bg_est_obj = Bg_est_obj .* Bg_mask; % Background estimated object masked by background mask (Func return val.)

        curvature = curvature_central(phi);
        delta_phi = Dirac_global(phi,2);

        % eq.7 in chan-vese segmentation paper BY: Pascal Getreuer, delta_phi is H(phi(x))
        evolve_force = delta_phi .* (-prm.lambda1 * (prm.Featured_Img - Fg_est_obj).^2 + prm.lambda2*(prm.Featured_Img - Bg_est_obj).^2);

        reg_force = prm.length_term * curvature; % Adjusted heaviside
 
        dphi_dt = evolve_force./(max(abs(evolve_force(:)))+eps) + reg_force;

        delta_t = 0.8/(max(abs(dphi_dt(:)))+eps); % Step size using CFL

        prev_mask = phi >=0;
        
        phi = phi + phi_update_rate * delta_t * dphi_dt;

        phi = SussmanReinitLS(phi, 0.5); % Reinitialize LSF by Sussman reinitialization method
        phi = NeumannBoundCond(phi); % Check Neumann boundary condition
        
        curr_mask = phi >=0 ;
        
        if prm.display_intrvl > 0
            if mod(its, prm.display_intrvl) == 0
                fig = findobj(0,'tag','creaseg');
                ud = get(fig,'userdata');
                set(ud.txtInfo1,'string',sprintf('iteration: %d',its),'color',[1 1 0]);
                showCurveAndPhi(phi, prm.Img, prm.contour_color);
                drawnow;
            end
        end

        count = test_convergence(prev_mask, curr_mask, prm.convg_error, count);
        if count <= prm.convg_count
            its = its + 1;
        else
            stop = 1;
            %fprintf('Algorithm converged at iteration %d. \n', its);
            %figure(1);
            %title(sprintf('Algorithm converged at iteration %d. \n', its));
        end
    end

%     if stop == 0
%         fprintf('iteration limited violated (%d iteration(s)) - the algorithm does not converge!\n', its);
%         figure(1);
%         title(sprintf('iteration limited violated (%d iteration(s)) - the algorithm does not converge!\n', its));
%     end

    showCurveAndPhi(phi, prm.Img, prm.contour_color);
end

%% computeSDF - COMPUTESDF Create the signed distance function from the binary image.  inside >= 0, outside < 0
function [SDF] = computeSDF(bwI)
    phi = bwdist(bwI) - bwdist(1-bwI) + im2double(bwI)- 0.5;
    SDF = -phi;
end

%% Convergence Test
function count = test_convergence(prev_mask, curr_mask, convg_error, count)
    diff = prev_mask - curr_mask;
    n_diff = sum(abs(diff(:)));
    if n_diff < convg_error
        count = count + 1;
    else
        count = 0;
    end
end

%% Compute curvature    
function K = curvature_central(phi)                       
    [phi_x, phi_y] = gradient(phi);                                  
    normDu = sqrt(phi_x.^2 + phi_y.^2 + 1e-10); % plus a small possitive number to avoid division by zero.
    Nx = phi_x ./ normDu;                                       
    Ny = phi_y ./ normDu;
    Nxx = gradient(Nx);                              
    [~,Nyy] = gradient(Ny);                              
    K = Nxx + Nyy; % compute divergence
end

%% Check Neumann boundary condition
function g = NeumannBoundCond(f)
    [nrow, ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
end

%% SUSSMANREINITLS Reinitialize LSF by Sussman reinitialization method
function [D] = SussmanReinitLS(D,dt)
    %D  : level set function
    %dt : small timestep ~ 0.5
    
    shiftL = @(M) ([ M(:,2:size(M,2)) M(:,size(M,2)) ]);
    shiftR = @(M) ([ M(:,1) M(:,1:size(M,2)-1) ]);
    shiftD = @(M) (shiftR(M')');
    shiftU = @(M) (shiftL(M')');
    sussman_sign = @(D) (D ./ sqrt(D.^2 + 1));

    a = D - shiftR(D); % backward
    b = shiftL(D) - D; % forward
    c = D - shiftD(D); % backward
    d = shiftU(D) - D; % forward

    a_p = a;  a_n = a; % a+ and a-
    b_p = b;  b_n = b;
    c_p = c;  c_n = c;
    d_p = d;  d_n = d;

    a_p(a < 0) = 0;
    a_n(a > 0) = 0;
    b_p(b < 0) = 0;
    b_n(b > 0) = 0;
    c_p(c < 0) = 0;
    c_n(c > 0) = 0;
    d_p(d < 0) = 0;
    d_n(d > 0) = 0;

    dD = zeros(size(D));
    D_neg_ind = find(D < 0);
    D_pos_ind = find(D > 0);
    dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
    dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

    D = D - dt .* sussman_sign(D) .* dD;
end

%% showCurveAndPhi
% function showCurveAndPhi(phi, Img, contour_color)
% 	imshow(uint8(Img),[],'InitialMagnification',200);
%     hold on; 
% 	[c,~] = contour(phi,[0 0],contour_color,'Linewidth',2); hold off;
% 	test = isequal(size(c,2),0);
% 	while (test==false)
%         s = c(2,1);
%         if ( s == (size(c,2)-1) )
%             t = c;
%             hold on; plot(t(1,2:end)',t(2,2:end)',contour_color,'Linewidth',2);
%             test = true;
%         else
%             t = c(:,2:s+1);
%             hold on; plot(t(1,1:end)',t(2,1:end)',contour_color,'Linewidth',2);
%             c = c(:,s+2:end);
%         end
% 	end    
% end

function showCurveAndPhi(phi,ud,cl)
    fig = findobj(0,'tag','creaseg');
    ud = get(fig,'userdata');
	axes(get(ud.imageId,'parent'));
	delete(findobj(get(ud.imageId,'parent'),'type','line'));
	hold on; [c,h] = contour(phi,[0 0],cl{1},'Linewidth',2); hold off;
	delete(h);
    test = isequal(size(c,2),0);
	while (test==false)
        s = c(2,1);
        if ( s == (size(c,2)-1) )
            t = c;
            hold on; plot(t(1,2:end)',t(2,2:end)',cl{1},'Linewidth',3);
            test = true;
        else
            t = c(:,2:s+1);
            hold on; plot(t(1,1:end)',t(2,1:end)',cl{1},'Linewidth',3);
            c = c(:,s+2:end);
        end
    end
end




%% Generate the orthonormal 2d basis, compute K shifted 2d basis for the vectorized image
function [B] = Gen_2D_Vec_Basis_Fucn(Featured_Img, poly_degree, display)
    [Nr, Nc, ~] = size(Featured_Img);

    N = length(Featured_Img(:));     % Vectorized image

    B = zeros(N,(poly_degree+1).^2);

    [Basis_Func_row, ~] = Gen_1D_Vec_Basis_Fucn('revKrawtchouk', Nr, poly_degree); % revKrawtchouk legendre
    [Basis_Func_col, ~] = Gen_1D_Vec_Basis_Fucn('revKrawtchouk', Nc, poly_degree);

    ind = 0;
    for ii = 1 : poly_degree+1
        for jj = 1 : poly_degree+1
            ind = ind+1;
            row_basis = Basis_Func_row(:,ii);
            col_basis = Basis_Func_col(:,jj);
            outer_prod = row_basis * col_basis';  % same size as the image
            B(:,ind) = outer_prod(:);
        end
    end
    
    if (display == 1)
        figure(2);
        for i = 1 : poly_degree+1
            
            subplot(poly_degree+1, 2, 2*i-1);
            plot(Basis_Func_row(:,i));
            title(sprintf('1D Basis Func. [row#%d]',i));
            
            subplot(poly_degree+1, 2, 2*i);
            plot(Basis_Func_col(:,i));
            title(sprintf('1D Basis Func. [col#%d]',i));
        end
        
        figure(3)
        for i = 1 : (poly_degree+1)^2
            subplot(poly_degree+1, poly_degree+1, i)
            mesh(reshape(B(:,i),size(Featured_Img)))
            title(sprintf('2D Basis Func #%d',i));
        end

        figure(1)
    end
end

%% Gen_1D_Vec_Basis_Fucn; for 1D legendre see Kale and Vaswani
function [Basis_Func, orthonormal_Basis_Func] = Gen_1D_Vec_Basis_Fucn(Basis_Func, N, k)
    X = -1 : 2/(N-1) : 1;

    if strcmp(Basis_Func, 'legendre')
        p0 = ones(1,N);

        Basis_Func = zeros(N,k+1);
        orthonormal_Basis_Func = Basis_Func;
        Basis_Func(:,1) = p0';
        orthonormal_Basis_Func(:,1) = Basis_Func(:,1)/norm(Basis_Func(:,1));

        for ii = 2 : k+1
            Pn = 0;
            n = ii-1;   % degree
            for k = 0 : n
               Pn = Pn +  (nchoosek(n,k)^2)*(((X-1).^(n-k)).*(X+1).^k);
            end
            Basis_Func(:,ii) = Pn'/(2)^n;
            orthonormal_Basis_Func(:,ii) = Basis_Func(:,ii)/norm(Basis_Func(:,ii));
        end
        

    elseif strcmp(Basis_Func, 'revKrawtchouk')
        Basis_Func = zeros(N, k+1);
        orthonormal_Basis_Func = Basis_Func;
        
        % produce regressor
        Basis_Func = revKrawtchouk (k+1, 0.3, X, 2);

        for ii = 2 : k+1
            orthonormal_Basis_Func(:,ii) = Basis_Func(:,ii)/norm(Basis_Func(:,ii));
        end
        
    end
end

%% Krawtchouk
function X = revKrawtchouk (n, p, u, m)
    
    if ( n < 0 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'KRAWTCHOUK - Fatal error!\n' );
        fprintf ( 1, '  0 <= N is required.\n' );
        error ( 'KRAWTCHOUK - Fatal error!' );
    end

    if ( p <= 0.0 || 1.0 <= p )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'KRAWTCHOUK - Fatal error!\n' );
        fprintf ( 1, '  0 < P < 1 is required.\n' );
        error ( 'KRAWTCHOUK - Fatal error!' );
    end

    if ( m < 0 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'KRAWTCHOUK - Fatal error!\n' );
        fprintf ( 1, '  0 <= M is required.\n' );
        error ( 'KRAWTCHOUK - Fatal error!' );
    end
    
    %----------------------------------------------------------------------

    [~,h]=size(u);
    X = zeros(h, n+1);
    for k = 1:h
        X(k,1) = 1.0;
        X(k,2) = u(k) - (p * m);
    end

    for i = 1 : n - 1
        for k = 1:h
            X(k,i+2) = ( ( u(k) - ( i + p * ( m - 2 * i ) ) ) * X(k,i+1) ...
                - ( m - i + 1 ) * p * ( 1 - p ) * X(k,i)  ) ...
                / ( i + 1 );
        end
    end
end


%% EOF