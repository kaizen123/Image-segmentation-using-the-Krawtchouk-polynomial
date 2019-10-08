function [seg,evolved_phi,its] = creaseg_L2S(Img,init_mask,max_its ,length_term,poly_degree,l2_reg_term,color,thresh,display)
%      [bin_seg,convg_phi,its] = runL2S(Img,init_mask,max_iter,nu         ,m          ,lambda     ,color ,convg_thresh);

%CREASEG_SUVADIP: Version 1 of Smooth Basis Chan Vese 

%------------------------------------------------------------------------
% Legendre Basis Chan Vese: Faster implementation with vectorized basis
%
%
% Inputs: I           2D image
%         mask        Initialization (1 = foreground, 0 = bg)
%         max_its     Number of iterations to run segmentation for
%         length_term (optional) Weight of smoothing term
%                       higer = smoother.  default = 0.2
%         poly_degree Degree of the Legendre basis polynomial  
%         display     (optional) displays intermediate outputs
%                       default = true
%
% Outputs: seg        Final segmentation mask (1=fg, 0=bg)
%
% Description: The internal and external intensities are modeled as linear
% combination of legendre basis function, making it more adaptive to
% intensity inhomogeneity. Furthermore, the basis coefficients can be
% obtained by a close form solution, which makes it computationally
% feasible. This is the computationally faster version, where the basis
% vectors are vectorized.
% -------------------------------------------------------------------------
% IF YOU ARE USING THIS CODE, PLEASE CITE THE FOLLOWING PAPER:
% Mukherjee, S.; Acton, S., "Region based segmentation in presence of intensity 
% inhomogeneity using Legendre polynomials," Signal Processing Letters, IEEE , vol.PP, no.99, pp.1,1
% doi: 10.1109/LSP.2014.2346538

% phi = in inside>0 outside<0  bourder=0
phi                     = computeSDF(init_mask);
% phi1 = (phi-min(phi(:)))/(max(phi(:))-min(phi(:)))*255;
% class(phi)
% max(phi(:))
% min(phi(:))
% % figure(3),imshow(uint8(phi1));figure(1)
% figure(3),mesh(phi);figure(1)

n_poly                  = poly_degree ; 
[vect_Bases]            = LegendreBasis2D_vectorized(Img,n_poly);

% vect_Bases2 = vect_Bases;
% for i= 1+1:size(vect_Bases,1)-1
% for j= 1+1:size(vect_Bases,2)-1
%     k1 = vect_Bases(i-1:i+1,j-1:j+1);
%     vect_Bases2(i,j)=mean(k1(:));
% end
% end
% vect_Bases = vect_Bases2;

param.basis_vect        = vect_Bases;
param.init_phi          = phi;
param.Img               = im2graydouble(Img);

param.num_iter          = max_its;
param.convg_error       = thresh;
param.length_term       = length_term;
param.display_intrvl    = 30;

param.lambda1           = 1;
param.lambda2           = 1;
param.l2_reg_term       = l2_reg_term;

param.convg_count       = 20;
param.contour_color     = color;

[evolved_phi,its] = ChanVeseLegendre_vectorized(param);

seg = evolved_phi >= 0;
evolved_phi = -evolved_phi;     % Return the negative following Creaseg's convention

    % by kambiz rahbar
    save_resImage('L2S');

end

%--------------------- Function Heart ----------------------------------
function [phi,p1_in,p2_out] = ChanVeseLegendre_vectorized(opt)
%Segmentation by Legendre Level Set (L2S)

% Dirac_global is defrencial of Heaviside function
Dirac_global            = @(x,e) ((e/pi)./(e^2.+ x.^2));
Heaviside               = @(y,e) (0.5 * (1+ (2/pi) * atan(y/e) ) );

its         = 0;
max_iter    = opt.num_iter;
u           = opt.Img;
phi         = opt.init_phi;
convg_err   = opt.convg_error;
reg_term    = opt.length_term; %= nu 
l1          = opt.lambda1;   % =1
l2          = opt.lambda2;   % =1

display     = opt.display_intrvl;
color       = opt.contour_color;
count_lim   = opt.convg_count; % =20

% B is 2D Legendre
B           = opt.basis_vect;             % each col is the basis vector of length = vec(Img)
lambda_l2   = opt.l2_reg_term;            % small constant for L2 regularization to make matrix invertible

stop = 0;
count = 0;

II = eye(size(B,2),size(B,2));
while (its < max_iter && stop == 0)
    
    % h_phi is MASK
    h_phi = Heaviside(phi,2);
%     figure(4),mesh(h_phi);figure(1)
   
%     disp('size h_phi:')
%     size(h_phi)
%     disp('size phi:')
%     size(phi)
% max(h_phi(:))
% min(h_phi(:))

    % inside_mask is object MASK   and      outside_mask is background MASK
    inside_mask = h_phi;
    outside_mask = 1-h_phi;
    inside_indicator = inside_mask(:);
    outside_indicator = outside_mask(:);
%  disp('size U')
%  size(u)
%  disp('size inside mask')
%  size(inside_mask)
 
    % u_in & u_out are  MASKed image
    u_in = u.*inside_mask;
    u_out = u.*outside_mask;
    u_in = u_in(:);
    u_out = u_out(:);
    
    % A1 is  Legendre
    A1 = B';            % (each row contains a basis vector)
%     disp('size A1')
%     size(A1)
%     size(repmat(inside_indicator',size(A1,1),1))
   

    % A2 & B2 are  masked Legendre 
    A2 = A1.*(repmat(inside_indicator',size(A1,1),1));   % A1, with each row multiplied with hphi
    B2 = A1.*(repmat(outside_indicator',size(A1,1),1));   % A1, with each row multiplied with hphi
    
    % lambda_l2 is lambda in demo  =  0.1
%     disp('size A1*A2''')
%     size(II)

% use LSE to calculate c1_vec & c2_vec
    c1_vec = pinv(A1*A2' + lambda_l2*II) * A1*u_in;
    c2_vec = pinv(A1*B2' + lambda_l2*II) * A1*u_out;
    
%     c1_vec = (A1*A2' + lambda_l2*II)\(A1*u_in);
%     c2_vec = (A1*B2' + lambda_l2*II)\(A1*u_out);
%   disp('c1_vec')
%   size(c1_vec)

% calculate A'*p(x)  in   equation   #3
% p1_vec & p2_vec are estimation of object & background
    p1_vec = B*c1_vec;
    p2_vec = B*c2_vec;

    p1 = reshape(p1_vec,size(u));
    p2 = reshape(p2_vec,size(u));

% figure(6)
% subplot(221)
% mesh(p1)
% 
% subplot(222)
% mesh(p2)
% figure(1)

% p1_in: calculate A'*p(x) * h_phi  in   equation   #3

    p1_in = p1.*h_phi;
    p2_out = p2.*(1-h_phi);
    
    % calc diff2(phi)  phi is weighted init_mask  ,curvature = diff(phi,2,1) + diff(phi,2,2);
    curvature   = curvature_central(phi);
    % Dirac is diff of Heaviside
    delta_phi   = Dirac_global(phi,2);
%     mesh(delta_phi)
%     pause()

% equation 7 in chan-vese segmentation paper BY:Pascal Getreuer
% delta_phi is H(phi(x))
    evolve_force = delta_phi.*(-l1*(u-p1).^2 + l2*(u-p2).^2);
%     disp('size eval')
%     size(evolve_force)
%     evolve_force = h_phi.*(-l1*(u-p1).^2 + l2*(u-p2).^2);
    
% reg_force   is nu * diff2(phi)   taadil shodeye Heaviside
    reg_force    = reg_term*curvature;
%      max(abs(evolve_force(:)))
    dphi_dt = evolve_force./(max(abs(evolve_force(:)))+eps) + reg_force;
%    disp('---')
%    max(evolve_force(:))
%    max(abs(reg_force(:)))
    delta_t = .8/(max(abs(dphi_dt(:)))+eps);          % Step size using CFL
    
    
    prev_mask = phi >=0;
    
    phi = phi + delta_t*dphi_dt;
    
    %SUSSMANREINITLS Reinitialize LSF by Sussman reinitialization method
   phi = SussmanReinitLS(phi,0.5);
    phi = NeumannBoundCond(phi);
    if display > 0
%         if mod(its,display) == 0
%             showCurveAndPhi(phi,u,color);
%             drawnow;
%         end
    fig = findobj(0,'tag','creaseg');
    ud = get(fig,'userdata');
       if (display>0)
            if ( mod(its,15)==0 )            
                set(ud.txtInfo1,'string',sprintf('iteration: %d',its),'color',[1 1 0]);
                showCurveAndPhi(phi,ud,color);
                drawnow;
            end
        else
            if ( mod(its,10)==0 )            
                set(ud.txtInfo1,'string',sprintf('iteration: %d',its),'color',[1 1 0]);
                drawnow;
            end
        end


    end
    
    curr_mask = phi >=0 ;
    
    count = convergence(prev_mask,curr_mask,convg_err,count);
    % count how many succesive times we have attained convergence, reduce local minima
    if count <= count_lim
        its = its + 1;
    else
        stop = 1;
        fprintf('Algorithm converged, iteration=%d \n',its);
    end
    
end

% if stop == 0
%     fprintf('End of iteration %d',its);
% end

showCurveAndPhi(phi,u,color);
end


%--------------------- Auxilary Functions

function [SDF] = computeSDF(bwI)
%COMPUTESDF Create the signed distance function from the binary image
% inside >= 0, outside < 0

phi = bwdist(bwI)-bwdist(1-bwI)+im2double(bwI)- 0.5;
SDF = -phi ;

end



% Convergence Test
function c = convergence(p_mask,n_mask,thresh,c)
    diff = p_mask - n_mask;
    n_diff = sum(abs(diff(:)));
    if n_diff < thresh
        c = c + 1;
    else
        c = 0;
    end
end


% Compute curvature    
function k = curvature_central(u)                       

    [ux,uy] = gradient(u);                                  
    normDu = sqrt(ux.^2+uy.^2+1e-10);	% the norm of the gradient plus a small possitive number 
                                        % to avoid division by zero in the following computation.
    Nx = ux./normDu;                                       
    Ny = uy./normDu;
    nxx = gradient(Nx);                              
    [~,nyy] = gradient(Ny);                              
    k = nxx+nyy;                        % compute divergence
%     figure(6)
%     max(k(:))
%     min(k(:))
% k1 = (k-min(k(:)))/(max(k(:))-min(k(:)))*255;
% 
%     imshow(uint8(k1))
%     figure(1)
% mesh(k)
% pause()
end


% Check boundary condition
function g = NeumannBoundCond(f)
    
    [nrow,ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
end





% Generate the orthonormal Legendre bases (vectorized, see Kale and Vaswani)
function [B] = LegendreBasis2D_vectorized(Img,k)

%LEGENDREBASIS compute K shifted legendre basis for the vectorized image

% k is n_poly
[Nr,Nc] = size(Img);
N = length(Img(:));     % Vectorized image

B = zeros(N,(k+1).^2);
[B_r,B_r_ortho] = legendre_1D(Nr,k);
[B_c,B_c_ortho] = legendre_1D(Nc,k);

% [B_r,B_r_ortho] = my_func(Nr,k);
% [B_c,B_c_ortho] = my_func(Nc,k);

% size(B_r)
% disp('size b_r_ortho')
% size(B_r_ortho)
% figure(3)
% subplot(221)
% plot(B_r_ortho(:,1))
% subplot(222)
% plot(B_r_ortho(:,2))
% figure(1)

ind = 0;
for ii = 1 : k+1
    for jj = 1 : k+1
        ind = ind+1;
        row_basis = B_r(:,ii);
        col_basis = B_c(:,jj);
        outer_prod = row_basis*col_basis';  % same size as the image
        B(:,ind) = outer_prod(:);
        
    end
end
% disp('size B:')
% size(B)
% figure(2)
% for i= 1:9
%     subplot(3,3,i)
%     mesh(reshape(B(:,i),size(Img)))
% end
% 
% figure(1)
end





function [B,orthonormal_B] = legendre_1D(N,k)

X = -1:2/(N-1):1;
% disp('size X')
% size(X)
% disp(' K')
% 
%  k
p0 = ones(1,N);

B = zeros(N,k+1);
orthonormal_B = B;
B(:,1) = p0';
orthonormal_B(:,1) = B(:,1)/norm(B(:,1));

for ii = 2 : k+1
    Pn = 0;
    n = ii-1;   % degree
    for k = 0 : n
       Pn = Pn +  (nchoosek(n,k)^2)*(((X-1).^(n-k)).*(X+1).^k);
    end
    B(:,ii) = Pn'/(2)^n;
    orthonormal_B(:,ii) = B(:,ii)/norm(B(:,ii));
end
% figure(2)
% subplot(221)
% plot(B(:,1))
% subplot(222)
% plot(B(:,2))
% figure(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,orthonormal_B] = my_func(N,k)

X = -1:2/(N-1):1;

B = zeros(N,k+1);
B(:,1) = ones(1,N);
B(:,2) = X;
B(:,3) = X;



for ii = 2 : k+1
    orthonormal_B(:,ii) = B(:,ii)/norm(B(:,ii));
end
% figure(2)
% subplot(221)
% plot(B(:,1))
% subplot(222)
% plot(B(:,2))
% figure(1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D] = SussmanReinitLS(D,dt)
%SUSSMANREINITLS Reinitialize LSF by Sussman reinitialization method
%D  : level set function
%dt : small timestep ~ 0.5
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

function shift = shiftD(M)
    shift = shiftR(M')';
end

function shift = shiftL(M)
    shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];
end

function shift = shiftR(M)
    shift = [ M(:,1) M(:,1:size(M,2)-1) ];
end

function shift = shiftU(M)
    shift = shiftL(M')';
end

function S = sussman_sign(D)
    S = D ./ sqrt(D.^2 + 1);    
end



%-- Displays the image with curve superimposed
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
    


% function showCurveAndPhi(phi,Img,cl) % ORIGINAL L2S
% 
% 	imshow(Img,[],'InitialMagnification',200);
%     hold on; 
% 	[c,h] = contour(phi,[0 0],cl,'Linewidth',3); hold off;
% 	test = isequal(size(c,2),0);
% 	while (test==false)
%         s = c(2,1);
%         if ( s == (size(c,2)-1) )
%             t = c;
%             hold on; plot(t(1,2:end)',t(2,2:end)',cl,'Linewidth',3);
%             test = true;
%         else
%             t = c(:,2:s+1);
%             hold on; plot(t(1,1:end)',t(2,1:end)',cl,'Linewidth',3);
%             c = c(:,s+2:end);
%         end
% 	end    
% end



function img = im2graydouble(img)    
    [dimy, dimx, c] = size(img);
    if(isfloat(img)) % image is a double
        if(c==3) 
            img = rgb2gray(uint8(img)); 
        end
    else           % image is a int
        if(c==3) 
            img = rgb2gray(img); 
        end
        img = double(img);
    end
end

