function [y,jaco] = map(x,varargin)

    % map(x) returns a nxd matrix which contains n rows of
    % points coordinates resulting from mapping points contained in x.
    %
    % map(x,'mapType',mapType) specify which map to be used: Henon(default),
    % Zaslavsk, Losi, Standard, Lorenz
    %
    % map(x,'param') specify the parameter in the systems
    %
    
    p = inputParser;
    addRequired(p,'x');
    defaultType = 'Hen';
    defaultParam = zeros(size(x,2),1);
    defaultJacob = 0;
    defaultXn = '';
    defaultYn = '';
    defaultZn = '';
    defaultInlineon = 1;
    defaulty1 = '';
    defaulty2 = '';
    defaulty3 = '';
    defaultsigma = 0;
    addParameter(p,'jacob',defaultJacob,@isnumeric)
    addParameter(p,'mapType',defaultType,@ischar);
    addParameter(p,'param',defaultParam,@isvector);
    addParameter(p,'xn',defaultXn,@ischar);
    addParameter(p,'yn',defaultYn,@ischar);
    addParameter(p,'zn',defaultZn,@ischar);
    addParameter(p,'inline',defaultInlineon,@isnumeric);
    addParameter(p,'y1',defaulty1);
    addParameter(p,'y2',defaulty2);
    addParameter(p,'y3',defaulty3);
    addParameter(p,'sigma',defaultsigma);
    parse(p,'x',varargin{:});
    mapType = p.Results.mapType;
    param = p.Results.param;
    jacob = p.Results.jacob;
    x_n = p.Results.xn;
    y_n = p.Results.yn;
    z_n = p.Results.zn;
    inlineon = p.Results.inline;
    y1 = p.Results.y1;
    y2 = p.Results.y2;
    y3 = p.Results.y3;
    sigma = p.Results.sigma;
    
    if inlineon == 0
        y_1 = y1;
        y_2 = y2;
        y_3 = y3;
    end
    jaco = 0;
    n = size(x,1);
    d = size(x,2);
    y = zeros(n,d);
    for i = 1:n
        % linear map with small non-linear perturbation
        if mapType == 'noL'
            y_1 = param(1)*x(i,1) + param(2)*x(i,2) + sigma*x(i,1).^2;
            y_2 = param(3)*x(i,1) + param(4)*x(i,2) + sigma*x(i,2).^2;
            y(i,:) = [y_1 y_2];
        % linear map
        elseif mapType == 'Lin'
            y_1 = param(1)*x(i,1) + param(2)*x(i,2);
            y_2 = param(3)*x(i,1) + param(4)*x(i,2);
            y(i,:) = [y_1 y_2];
            
        % Zaslavsk
        elseif mapType == 'Zas'
            y_2 = cos(2*pi*x(i,1)) + exp(-3)*x(i,2);
            y_1 = mod(x(i,1) + param(2) + param(1)*y_2,1);
            y(i,:) = [y_1 y_2]; 
            if i == n && jacob
                y_2 = [-2*pi*sin(2*pi*x(i,1)),exp(-3)];
                y_1 = mod(param(1)*[1-2*pi*sin(2*pi*x(i,1)),exp(-3)],1);
                jaco = [y_1;y_2];
            end

            
        % Losi
        elseif mapType == 'Los' 
            y_1 = 1 + x(i,2)-param(1)*abs(x(i,1));
            y_2 = param(2)*x(i,1);
            y(i,:) = [y_1 y_2]; 
            if i == n && jacob
                y_1 = [-sign(x(i,1))*param(1) 1];
                y_2 = [param(2) 0];
                jaco = [y_1;y_2];
            end
        
        % classic henon
        elseif mapType == 'Hen'
            y_1 = 1 + x(i,2)-  param(1)*x(i,1)^2;
            y_2 = param(2)*x(i,1);
            %y_1 = sign(y_1)*mod(abs(y_1),3);
            %y_2 = sign(y_2)*mod(abs(y_2),3);
            y(i,:) = [y_1 y_2]; 
            if i == n && jacob
                y_1 = [2*param(1)*x(i,1),1];
                y_2 = [param(2) 0];
                jaco = [y_1;y_2];
            end
        
        % inverse henon
        elseif mapType == 'iHe'
            y_1 = x(i,2)/param(2);
            y_2 = x(i,1) - 1 - param(1)*x(i,2)^2/param(2)^2;
            %y_1 = sign(y_1)*mod(abs(y_1),3);
            %y_2 = sign(y_2)*mod(abs(y_2),3);
            y(i,:) = [y_1 y_2]; 
            if i == n && jacob
                y_1 = [2*param(1)*x(i,1),1];
                y_2 = [param(2) 0];
                jaco = [y_1;y_2];
            end
        
        % standard map
        elseif mapType == 'Sta'
            y_1 = x(i,1) + param(1)*sin(x(i,2));
            y_2 = param(2)*(x(i,2)+y_1);
            y(i,:) = [y_1 y_2]; 
            if i == n && jacob
                y_1 = [1,param(1)*cos(x(i,2))];
                y_2 = param(2)*[1,1+param(1)*cos(x(i,2))];
                jaco = [y_1;y_2];
            end
        
        % lorenz
        elseif mapType == 'Lor'
            y_1 = param(1)*(x(i,2)-x(i,1));
            y_2 = param(2)*x(i,1) - x(i,1)*x(i,3) - x(i,2);
            y_3 = x(i,1)*x(i,2) - 8/3*x(i,3);
            y(i,:) = [y_1 y_2 y_3]; 
            if i == n && jacob
                y_1 = [-param(1),param(1),0];
                y_2 = [param(2)-x(i,3), -1, -x(i,1)];
                y_3 = [x(i,2), x(i,1),-8/3];
                jaco = [y_1;y_2;y_3];
            end
        
        elseif mapType == 'Opt'
            if d == 2
                if inlineon
                    y_1 = str2func(['@(x,y,a,b)',x_n]);
                    y_2 = str2func(['@(x,y,a,b)',y_n]);
                end
                y_11 = y_1(x(i,1),x(i,2),param(1),param(2));
                y_22 = y_2(x(i,1),x(i,2),param(1),param(2));
                y(i,:) = [y_11 y_22]; 
            elseif d == 3
                if length(param) == 2
                    if inlineon
                        y_1 = str2func(['@(x,y,z,a,b)',x_n]);
                        y_2 = str2func(['@(x,y,z,a,b)',y_n]);
                        y_3 = str2func(['@(x,y,z,a,b)',z_n]);
                    end
                    y_11 = y_1(x(i,1),x(i,2),x(i,3),param(1),param(2));
                    y_22 = y_2(x(i,1),x(i,2),x(i,3),param(1),param(2));
                    y_33 = y_3(x(i,1),x(i,2),x(i,3),param(1),param(2));
                    y(i,:) = [y_11 y_22 y_33];
                elseif length(param) == 3
                    if inlineon
                        y_1 = str2func(['@(x,y,z,a,b,c)',x_n]);
                        y_2 = str2func(['@(x,y,z,a,b,c)',y_n]);
                        y_3 = str2func(['@(x,y,z,a,b,c)',z_n]);
                    end
                    y_1 = y_1(x(i,1),x(i,2),x(i,3),param(1),param(2),param(3));
                    y_2 = y_2(x(i,1),x(i,2),x(i,3),param(1),param(2),param(3));
                    y_3 = y_3(x(i,1),x(i,2),x(i,3),param(1),param(2),param(3));
                    y(i,:) = [y_1 y_2 y_3];
                end
            end
        elseif mapType == 'Com'
            z = x(i,1) + x(i,2)*1i;
            zz = z-(z.^3-2*z+2)./(3*z.^2-2);
            y_1 = real(zz);
            y_2 = imag(zz);
            %y_1 = sign(y_1)*mod(abs(y_1),3);
            %y_2 = sign(y_2)*mod(abs(y_2),3);
            y(i,:) = [y_1 y_2];    
        end
    end
end