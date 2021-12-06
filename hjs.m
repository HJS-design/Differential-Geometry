classdef hjs < matlab.System
    % Untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    % Public, tunable properties
    properties
        n;
        coordinate;
        X;
        g;
        connection;
        Riemann_curvature_tensor;
        Riemann_curvature_metric_tensor;
        Ricci_tensor;
        Scalar_curvature;
        Gauss_curvature;
        
        alpha=1;
        scale=0.9;
    end
    
    
    
    methods
        function obj = hjs(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
        function []= geodesic(obj,interval,Init,interval_num)
            syms u(t)  T(t) [2,1] ;
            if nargin < 4
                ninterval=100;
            else
                ninterval=interval_num;
            end
            
            Y=sym('Y',[4,1]);
            coordinate=obj.coordinate;
            n=obj.n;
            X=obj.X;
            ss=obj.connection;
            sst=subs(ss,coordinate,u);
            dR=diff(u,1);
            T=dR;
            s1=diff(T,1)+squeeze(sum(repmat(reshape(T,1,n,1),[n,1,n]).* ...
                repmat(reshape(dR,1,1,n),[n,n,1]).*sst,[2,3]))==0;
            [V,S]=odeToVectorField((s1));
            M=matlabFunction(V,'vars',{'t','Y'});
            syms u1 u2  Du1 Du2;
            yInit=matlabFunction(subs([u1;u2;Du1;Du2],S,Y));
            
            ySol = ode45(M,interval,yInit(Init(1),Init(2),Init(3),Init(4)));
            tValues = linspace(interval(1),interval(2),ninterval);
            y=yInit(1 ,2 ,3, 4);
            yValues = deval(ySol,tValues,y(1:2));
            y=matlabFunction(X);
            tt=y(yValues(1,:),yValues(2,:));
            
            if(size(tt,1)==2)
                plot(tt(1,:),tt(2,:));
            else
                plot3(tt(1,:),tt(2,:),tt(3,:),"LineWidth",2);
            end
            
            
            
        end
        function [yValues]= parallel_transport(obj,u,interval,Init,interval_num)
            syms  T(t) [2,1] ;
            if nargin < 5
                ninterval=100;
            else
                ninterval=interval_num;
            end
            
            n=obj.n;
            coordinate=obj.coordinate;
            X=obj.X;
            [V,S]=odeToVectorField(diff(T,1)+squeeze(sum(repmat(reshape(T,1,n,1),[n,1,n]).*repmat(reshape(diff(u,1),1,1,n),[n,n,1]).*subs(obj.connection,coordinate,u),[2,3]))==0);
            M=matlabFunction(V,'vars',{'t','Y'});
            syms T1 T2;
            Y=sym('Y',[2,1]);
            yInit=matlabFunction(subs([T1;T2],S,Y));
            
            ySol = ode45(M,interval,yInit(Init(1), Init(2)));
            
            tValues = linspace(interval(1),interval(2),ninterval);
            yValues = deval(ySol,tValues,yInit(1 ,2));
           
            syms tt;
            temp(t,tt)=subs(X,coordinate,u)+repmat(tt,[size(X,1),1]);
            t2Xf=matlabFunction(temp,'vars',{'t','tt'});
            identity=matlabFunction(tt,'vars',{'tt'});
            Xvalue=t2Xf(tValues,tValues)-identity(tValues);
            temp(t,tt)=subs([diff(X,coordinate(1));diff(X,coordinate(2))],coordinate,u)+repmat(tt,[size(X,1)*2,1]);
            dXf=matlabFunction(temp,'vars',{'t','tt'});
            rdx=squeeze(sum(reshape(dXf(tValues,tValues)-identity(tValues),[size(X,1),2,ninterval]).*repmat(reshape(yValues,[1,2,ninterval]),[size(X,1),1,1]),2));
           
            if(size(X,1)==2)
                quiver(Xvalue(1,:),Xvalue(2,:),rdx(1,:),rdx(2,:),obj.scale,"LineWidth",2,'ShowArrowHead',0);
            else
                quiver3(Xvalue(1,:),Xvalue(2,:),Xvalue(3,:),rdx(1,:),rdx(2,:),rdx(3,:),obj.scale,"LineWidth",2,'ShowArrowHead',0);
            end
            hold on;
            if(size(X,1)==2)
                plot(Xvalue(1,:),Xvalue(2,:));
            else
                plot3(Xvalue(1,:),Xvalue(2,:),Xvalue(3,:),"LineWidth",2);
            end
            hold off;
          
        end
        function []= drawmesh(obj,uinterval,vinterval,draw_curvature,interval_num)
            if(size(obj.X,1)==2)
                return;
            end
            if nargin < 4
                draw_curvature = 0;
            end
            if nargin < 5
                interval_num =   30;
            end
            xvalue=linspace(uinterval(1),uinterval(2),interval_num);
            yvalue=linspace(vinterval(1),vinterval(2),interval_num);
            [x,y]=ndgrid(xvalue,yvalue);
            Xf=matlabFunction(obj.X,'vars',{'u1','u2'});
            DT = delaunayTriangulation(x(:),y(:));
            [C,IA,IC]=uniquetol(Xf(DT.Points(:,1)',DT.Points(:,2)')','ByRows',true);
            TR=triangulation(IC(DT.ConnectivityList),C);
            VV=vertexNormal(TR);
            if(draw_curvature)
                curvatureF=matlabFunction(obj.Scalar_curvature,'vars',{'u1','u2'});
                color=curvatureF(DT.Points(:,1),DT.Points(:,2));
                if(numel(color)==1)
                    color=repmat(color,[size(DT.Points,1),1]);
                    
                end
                patch('Faces',IC(DT.ConnectivityList),'Vertices',C,'FaceVertexCData',color(IA),'FaceColor','interp','EdgeColor','none','VertexNormals',VV,'FaceLighting','gouraud','BackFaceLighting','unlit',"FaceAlpha",obj.alpha);
                colormap parula;
                view(3);
             
                colorbar;
            else
                patch('Faces',IC(DT.ConnectivityList),'Vertices',C,'FaceVertexCData',C(:,1),'FaceColor','interp','EdgeColor','none','VertexNormals',VV,'FaceLighting','gouraud','BackFaceLighting','unlit',"FaceAlpha",obj.alpha);
                camlight;
                view(3);
          
            end
            
            
        end
        function []= curvature_tensor(obj)
            n=obj.n;
            connection=obj.connection;
            g=obj.g;
            coordinate=obj.coordinate;
            t1=cat(4,diff(connection,coordinate(1)),diff(connection,coordinate(2)));
            t2=squeeze(sum(repmat(reshape(connection,1,n,n,n),[n,1,1,1,n]).*repmat(reshape(connection,n,n,1,1,n),[1,1,n,n,1]),2));
            t3=simplify(permute(t1,[1,2,4,3])-t1+permute(t2,[1,2,4,3])-t2);
            obj.Riemann_curvature_tensor=t3;
            
            t1=simplify(squeeze(sum(repmat(reshape(g,n,n),[1,1,n,n,n]).*repmat(reshape(t3,1,n,n,n,n),[n,1,1,1,1]),2)));
            obj.Riemann_curvature_metric_tensor=t1;
            
            [x1,~,x3,~]=ndgrid(1:n,1:n,1:n,1:n);
            t3(x1~=x3)=0;
            
            obj.Ricci_tensor=simplify(squeeze(sum(t3,[1,3])));
            
            ig=inv(g);
            
            obj.Scalar_curvature=simplify(sum(ig.*obj.Ricci_tensor,"all"));
            obj.Gauss_curvature=obj.Scalar_curvature/2;
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            
            
            coordinate=obj.coordinate;
            n=obj.n;
            if(isempty(obj.g))
                X=obj.X;
                df=[diff(X,coordinate(1)),diff(X,coordinate(2))];
                obj.g=simplify((df.')*df);
            end
            g=obj.g;
            clear df;
            ig=inv(g);
            dg=cat(3,diff(g,coordinate(1)),diff(g,coordinate(2)));
            s1=squeeze(sum(repmat(reshape(ig,n,1,n),[1,n,1,n]).*repmat(reshape(dg,1,n,n,n),[n,1,1,1]),3));
            s2=squeeze(sum(repmat(reshape(ig,n,1,1,n),[1,n,n,1]).*repmat(reshape(dg,1,n,n,n),[n,1,1,1]),4));
            obj.connection=(permute(s1,[1,3,2])+s1-s2)/2;
            
            % Perform one-time calculations, such as computing constants
        end
        
        function []= stepImpl(obj)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            
        end
        
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end