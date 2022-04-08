function []=Testing_SA()

CB=4; % cantidad de Instancias

NA=6; % cantidad de algoritmos

Co=30; % Cantidad de corridas

Par=6; % Cantidad de valores de Alpha

Alpha=[0.5 0.6 0.7 0.8 0.9 0.99];

T=[300,800];

CI=[5000,10000];

%Corridas

for j=1:CB % Para cada una de las instancias
    
    switch j
        case 1
            disp('5P30A')
           load('/home/bolivar/Escritorio/bolivar/Escritorio/SA/5p30a.mat')
        case 2
             disp('3P60A')
            load('/home/bolivar/Escritorio/bolivar/Escritorio/SA/3p60a.mat')
        case 3
             disp('2P120A')
            load('/home/bolivar/Escritorio/bolivar/Escritorio/SA/2p120a.mat')
        case 4
             disp('3P90A')
            load('/home/bolivar/Escritorio/bolivar/Escritorio/SA/3p90a.mat')
    end
        
for a=1:NA % Para cada uno de los algoritmos
switch  a
    
    case 1
       disp('SA1')
          
       Data=cat(3,zeros(Co,Par),zeros(Co,Par));
       
       Stat=zeros(2,4);
     
 for t=1:2      % Para los valores de temperatura y cantidad de itereaciones
         switch t
             case 1
                  aux1=100;
             case 2
                  aux2=100;
         end
            
         tic   
         for p=1:Par   % Para la cantidad de valores de Alpha
             
           for i=1:Co % Para la acantidad de corridas
                                               
             [FO,Data(i,p,t),BestSchedule] = SA01(DataProject,RRHH,TA,CP,CI(t),T(t),Alpha(p));
                                   
             switch t 
                 
                 case 1
                     
                    if  Data(i,p,t) < aux1 
                 
                         aux1= Data(i,p,t);
                 
                         Best1=BestSchedule;
                 
                         FB1=FO;                                      
                    end

                 case 2
                     
                     if  Data(i,p,t) < aux2  
                 
                         aux2= Data(i,p,t);
                 
                         Best2=BestSchedule;
                 
                         FB2=FO;                                      
                     end
                     
             end                          
           end
           
         end
          toc
          
       switch t
           case 1
            DAT=Data(:,:,t);
            Stat(1,1)=min(DAT(1:numel(DAT)));  Stat(1,2)=max(DAT(1:numel(DAT))); 
            Stat(1,3)=mean(DAT(1:numel(DAT))); Stat(1,4)=std(DAT(1:numel(DAT)));
           case 2
            DAT=Data(:,:,t);   
            Stat(2,1)=min(DAT(1:numel(DAT)));  Stat(2,2)=max(Data(1:numel(DAT))); 
            Stat(2,3)=mean(DAT(1:numel(DAT))); Stat(2,4)=std(Data(1:numel(DAT)));
       end
                    
 end   
       FB=[FB1,FB2]; 
       Best=struct('Best1',Best1,'Best2',Best2);
    switch j
      
       case 1
          save('SA1_5P30A','Data','FB','Best','Stat') 
      case 2
          save('SA1_3P60ASA1','Data','FB','Best','Stat') 
      case 3
          save('SA1_2P120ASA1','Data','FB','Best','Stat') 
      case 4
          save('SA1_3P90SA1','Data','FB','Best','Stat') 
    end
        
% Algoritmo SA 2    

    case 2
        disp('SA2')
          
         Data=cat(3,zeros(Co,Par),zeros(Co,Par));
       
       Stat=zeros(2,4);
     
 for t=1:2      % Para los valores de temperatura y cantidad de itereaciones
         switch t
             case 1
                 aux1=100;
             case 2
                  aux2=100;
         end
            
         tic   
         for p=1:Par   % Para la cantidad de valores de Alpha
             
           for i=1:Co  % Para la acantidad de corridas
                                               
             [FO,Data(i,p,t),BestSchedule] = SA02(DataProject,RRHH,TA,CP,W,CI(t),T(t),Alpha(p));
                                    
             switch t 
                 
                 case 1
                     
                    if  Data(i,p,t) < aux1 
                 
                         aux1= Data(i,p,t);
                 
                         Best1=BestSchedule;
                 
                         FB1=FO;                                      
                    end

                 case 2
                     
                     if  Data(i,p,t) < aux2  
                 
                         aux2= Data(i,p,t);
                 
                         Best2=BestSchedule;
                 
                         FB2=FO;                                      
                     end
                     
             end                          
           end
           
         end
          toc
          
       switch t
          case 1
            case 1
            DAT=Data(:,:,t);
            Stat(1,1)=min(DAT(1:numel(DAT)));  Stat(1,2)=max(DAT(1:numel(DAT))); 
            Stat(1,3)=mean(DAT(1:numel(DAT))); Stat(1,4)=std(DAT(1:numel(DAT)));
           case 2
            DAT=Data(:,:,t);   
            Stat(2,1)=min(DAT(1:numel(DAT)));  Stat(2,2)=max(Data(1:numel(DAT))); 
            Stat(2,3)=mean(DAT(1:numel(DAT))); Stat(2,4)=std(Data(1:numel(DAT)));
       end
                    
 end   
       FB=[FB1,FB2]; 
       Best=struct('Best1',Best1,'Best2',Best2);
   switch j
      case 1
          save('5P30ASA2','Data','FB','Best','Stat') 
      case 2
          save('3P60ASA2','Data','FB','Best','Stat') 
      case 3
          save('2P120ASA2','Data','FB','Best','Stat') 
      case 4
           save('3P90SA2','Data','FB','Best','Stat') 
   end
   
% Algoritmo SA3

    case 3
       disp('SA3')
       
       Data=cat(3,zeros(Co,Par),zeros(Co,Par));
       
       Stat=zeros(2,4);
     
 for t=1:2      % Para los valores de temperatura y cantidad de itereaciones
         switch t
             case 1
                 aux1=100;
             case 2
                  aux2=100;
         end
            
         tic   
         for p=1:Par   % Para la cantidad de valores de Alpha
             
           for i=1:Co  % Para la acantidad de corridas
                                               
             [FO,Data(i,p,t),BestSchedule] =  SA03(DataProject,RRHH,TA,CP,W,CI(t),T(t),Alpha(p));
                                    
             switch t 
                 
                 case 1
                     
                    if  Data(i,p,t) < aux1 
                 
                         aux1= Data(i,p,t);
                 
                         Best1=BestSchedule;
                 
                         FB1=FO;                                      
                    end

                 case 2
                     
                     if  Data(i,p,t) < aux2  
                 
                         aux2= Data(i,p,t);
                 
                         Best2=BestSchedule;
                 
                         FB2=FO;                                      
                     end
                     
             end                          
           end
           
         end
          toc
          
       switch t
           case 1
            DAT=Data(:,:,t);
            Stat(1,1)=min(DAT(1:numel(DAT)));  Stat(1,2)=max(DAT(1:numel(DAT))); 
            Stat(1,3)=mean(DAT(1:numel(DAT))); Stat(1,4)=std(DAT(1:numel(DAT)));
           case 2
            DAT=Data(:,:,t);   
            Stat(2,1)=min(DAT(1:numel(DAT)));  Stat(2,2)=max(Data(1:numel(DAT))); 
            Stat(2,3)=mean(DAT(1:numel(DAT))); Stat(2,4)=std(Data(1:numel(DAT)));
       end
                    
 end   
       FB=[FB1,FB2]; 
       Best=struct('Best1',Best1,'Best2',Best2); 
                         
   switch j
      case 1
          save('5P30ASA3','Data','FB','Best','Stat') 
      case 2
          save('3P60ASA3','Data','FB','Best','Stat') 
      case 3
          save('2P120ASA3','Data','FB','Best','Stat') 
      case 4
           save('3P90SA3','Data','FB','Best','Stat') 
   end
    
 % Algoritmo SA4
 
    case 4
         disp('SA4')
          
         Data=cat(3,zeros(Co,Par),zeros(Co,Par));
       
       Stat=zeros(2,4);
     
        for t=1:2      % Para los valores de temperatura y cantidad de itereaciones
         switch t
             case 1
                 aux1=100;
             case 2
                  aux2=100;
         end
            
         tic   
         for p=1:Par   % Para la cantidad de valores de Alpha
             
           for i=1:Co  % Para la acantidad de corridas
                                               
             [FO,Data(i,p,t),BestSchedule] =  SA04(DataProject,RRHH,TA,CP,W,CI(t),T(t),Alpha(p));
                                    
             switch t 
                 
                 case 1
                     
                    if  Data(i,p,t) < aux1 
                 
                         aux1= Data(i,p,t);
                 
                         Best1=BestSchedule;
                 
                         FB1=FO;                                      
                    end

                 case 2
                     
                     if  Data(i,p,t) < aux2  
                 
                         aux2= Data(i,p,t);
                 
                         Best2=BestSchedule;
                 
                         FB2=FO;                                      
                     end
                     
             end                          
           end
           
         end
          toc
          
       switch t
           case 1
            DAT=Data(:,:,t);
            Stat(1,1)=min(DAT(1:numel(DAT)));  Stat(1,2)=max(DAT(1:numel(DAT))); 
            Stat(1,3)=mean(DAT(1:numel(DAT))); Stat(1,4)=std(DAT(1:numel(DAT)));
           case 2
            DAT=Data(:,:,t);   
            Stat(2,1)=min(DAT(1:numel(DAT)));  Stat(2,2)=max(Data(1:numel(DAT))); 
            Stat(2,3)=mean(DAT(1:numel(DAT))); Stat(2,4)=std(Data(1:numel(DAT)));
       end
                    
        end   
       
       FB=[FB1,FB2]; 
       Best=struct('Best1',Best1,'Best2',Best2);
       
   switch j
      case 1
          save('5P30ASA4','Data','FB','Best','Stat') 
      case 2
          save('3P60ASA4','Data','FB','Best','Stat') 
      case 3
          save('2P120ASA4','Data','FB','Best','Stat') 
      case 4
           save('3P90SA4','Data','FB','Best','Stat') 
   end
 % Algoritmo SA5
    case 5
         disp('SA5')
          
         Data=cat(3,zeros(Co,Par),zeros(Co,Par));
       
       Stat=zeros(2,4);
     
        for t=1:2      % Para los valores de temperatura y cantidad de itereaciones
         switch t
             case 1
                 aux1=100;
             case 2
                  aux2=100;
         end
            
         tic   
         for p=1:Par   % Para la cantidad de valores de Alpha
             
           for i=1:Co  % Para la acantidad de corridas
                                               
             [FO, Data(i,p,t),BestSchedule] =  SA05(DataProject,RRHH,TA,CP,W,CI(t),T(t),Alpha(p));;
                                    
             switch t 
                 
                 case 1
                     
                    if  Data(i,p,t) < aux1 
                 
                         aux1= Data(i,p,t);
                 
                         Best1=BestSchedule;
                 
                         FB1=FO;                                      
                    end

                 case 2
                     
                     if  Data(i,p,t) < aux2  
                 
                         aux2= Data(i,p,t);
                 
                         Best2=BestSchedule;
                 
                         FB2=FO;                                      
                     end
                     
             end                          
           end
           
         end
          toc
          
       switch t
           case 1
            DAT=Data(:,:,t);
            Stat(1,1)=min(DAT(1:numel(DAT)));  Stat(1,2)=max(DAT(1:numel(DAT))); 
            Stat(1,3)=mean(DAT(1:numel(DAT))); Stat(1,4)=std(DAT(1:numel(DAT)));
           case 2
            DAT=Data(:,:,t);   
            Stat(2,1)=min(DAT(1:numel(DAT)));  Stat(2,2)=max(Data(1:numel(DAT))); 
            Stat(2,3)=mean(DAT(1:numel(DAT))); Stat(2,4)=std(Data(1:numel(DAT)));
       end
                    
       end   
       FB=[FB1,FB2]; 
       Best=struct('Best1',Best1,'Best2',Best2);
        
   switch j
      case 1
          save('5P30ASA5','Data','FB','Best') 
      case 2
          save('3P60ASA5','Data','FB','Best') 
      case 3
          save('2P120ASA5','Data','FB','Best') 
      case 4
           save('3P90SA5','Data','FB','Best') 
   end
                
end % Para el caso del algoritmo a
end % Para cada una de los algoritmos
end % Para cada una de las instancias
disp('EL EXPERIMENTO TERMINO EXITOSAMENTE')
end % FUNCTION
