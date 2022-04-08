function [FO,M,BestSchedule] = SA02(DataProject,RRHH,TA,CP,CI,To,Alpha)
% Simulate Annealing function 
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
% Facultad III. Universidad de las Ciencias Informáticas
% Descripción:
% Entradas
% DataProject
%
% CI -> Cantidad de Inteciones
% To -> 
% Salidas
% INICIALIZACION
%tic

[SolInit,FVo,NP,CA] = SSGSpdRNew(DataProject,RRHH,TA,CP);

Best=FVo;

i=0; 

M=1000;

T=To;

FO=zeros(1,CI);

 while i < CI 
    
      i=i+1;
     
      [Schedule]=NeighborhoodR(DataProject,SolInit,NP,CA);
     
      [Schedule,FV]=SSGSSA(Schedule,DataProject,RRHH,TA,CP,NP,CA);
    
     Deltha=FV-FVo;
     
     if Deltha>0
         
     Aceptance =exp(-(Deltha)/T);
     
     end
    
     Prob=rand;
     
     if Deltha <= 0 || Aceptance > Prob
        
         SolInit=Schedule;
        
         Best=FV;
         
         if  Best < M
             
             M=Best;
             
             BestSchedule=SolInit;
             
         end
         
     end
     
      FVo=Best;
     
      FO(i)=Best;
                      
     T=To/(i*Alpha);
          
 end
%toc
end

function [Schedule,FO,NP,CA] = SSGSpdRNew(DataProject,RRHH,TA,CP)
% SSGSpd Summary of this function [Schedule,FO,CA,NP] = SSGSpdRNew(DataProject,RRHH,TA,CP)
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
% Facultad III. Universidad de las Ciencias Informáticas
% Descripción: Heuristica para la generación de soluciones iniciales factibles basada
% en Serial Scheme Generation Sechdeule (SSGS); con enfoque Multi-Proyecto, 
% la prioridad entre los proyecto cambia en el tiempo, la prioridad entre las actividades 
% está dada por la regla de prioridad MTS-modificada  y prioridad para la
% asignación de los recursos.
% Diseñada para  n proyectos n=2,3,... con la misma cantidad de actividades
% con tiempos de arribo Start=[T1,T2,T3,T4,...,Tn]
% los proyectos son tomados con dos actividades ficticias
% ENTRADA
%  DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Información sobre los recursos 
%  TA -> Tiempos de arribo de los i proyectos  TA=[TA1,TA2,...,TAn]
%  CP -> Rutas criticas de los i proyectos     CP=[CP1,CP2,...,CPn]
% SALIDAS
%  Schedule -> Solución  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  TI -> Tiempos de Inicio
%  TF -> Tiempos de Finalización 
%  RA -> Recurso Asignado                                         
%  W  -> Prioridad (MTS)

% DEFINICIÓN DE VARIABLES
%tic

[CA,H,NP]=size(DataProject); % Datos de los Proyectos (CA -> cantidad de actividades, NP-> cantidad de Proyectos)

Schedule=cat(3,zeros(NP,CA),zeros(NP,CA),zeros(NP,CA),rand(NP,CA)); % Estructura de la solución

PS=zeros(1,NP); % Variable que guarda los proyectos que han comenzado PS(i)={0,1} 

PF=PS;      % Variable que guarda los proyectos que han finalizado PF(i)={0,1} 

APP=zeros(NP,CA);  % Variable que guarda las actividades que se han planificado APP(i,j)={0,1} 

% INICIALIZACIÓN

PS(1)=1; %  Indica que el proyecto 1 comienza procesandose

APP(1,1)=1; %  Se planifica la primera actividad del primer proyecto

ARP=APP; % Conjuntos de actividades que pueden ser planificadas se inicializa con  la actividad 1 del proyecto 1  
 
Stop=0; %  Se inicializa en cero la condición de parada

ind=2:NP; % indice de los proyectos 

ind2=1:NP;

v=1; % Contador de los proyectos que inician

RDE=[];

while Stop < NP % Mientras la cantidad de proyectos finalizado sea menor que la actidad de proyectos 
    
if v~=NP 

    MK=max(max(Schedule(:,2:end,2)));
            
for k=nonzeros(ind)' % Para cada uno de los Proyectos
    
    if TA(k)<=MK && TA(k)~=0 % Si el mayor si el Makespan parcial MK es mayor o igual a la fecha de arribo deque la fe
          
       PS(k)=1;
       
       Schedule(k,1,1)=TA(k); % se inician los proyectos que aun no han comenzado
       
       Schedule(k,1,2)=TA(k);
              
       v=v+1;
       
       ind(k-1)=0;
       
       ARP(k,1)=1;
       
       APP(k,1)=1;
                    
     end
       
end
    
end
       
    PR=find(PS==1);% Proyectos que han comenzado 
    
    CPS=length(PR); % Cantidad de proyectos iniciados
   
     if CPS==1 % Si la cantidad de proyectos iniciados es 1 
         
         OrdenP=PR; % Orden es PR
         
    else % en caso contrario
        
         OrdenP=zeros(1,length(PR));
         
         MKP=OrdenP;
         
                
         for l=PR
            
           MKP(l) = max(Schedule(l,:,2));
                                                            
         end
         
         aux=1-(nonzeros(MKP)').^(-1);
              
         ORD=sortrows([aux;PR]',1);
         
         OrdenP=ORD(:,2)';
         
    end

    for i=OrdenP % Para proyectos Iniciado segun OrdenP
                                                                 
         [Schedule,APP,ARP,RDO] =SFMtsRa(Schedule,DataProject,RRHH,i,nonzeros(ARP(i,:))',ARP,APP,RDE);                                                          
               
         RDE=RDO;             
  
    end
        
  for r=nonzeros(ind2)'
    
    if Schedule(r,CA,2)~=0 && sum(APP(r,:))==CA
        
        PF(r)=1;
        
        PS(r)=0;   
        
        ind2(r)=0;
        
        ARP(r,:)=0;
             
    end
    
  end
  
  Stop=sum(PF);
  
end  
FO=mean(((Schedule(:,CA,2)'-(CP+TA))));  %mean((Schedule(:,CA,2)'./CP)-1)); %mean(Schedule(:,CA,2));
%toc
end

function [Schedule,APP,ARP,RDO] = SFMtsRa(Schedule,DataProject,RRHH,i,ASched,ARP,APP,RDO)
%  SchedulingFeasible Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
%  Facultad III. Universidad de las Ciencias Informáticas
%  Planifica las activitades de ActReadyPlan sin violar las restricciones de
% Precedencia y de Recursos.
% ERNTRADAS
% Schedule -> solución 
% DataProyect -> Datos de los Proyectos 
%ActReadyPlan-> Conjunto de actividades que están lista para ser planificadas 
% RRHH-> Datos sobre los recursos humanos del proyecto
% APP(i) -> Conjunto de Actividades planificadas por Proyecto.
% DataProject -> Informaci�n sobre los proyectos.
% SALIDAS
% Schedule-> solución 
% APP(i) -> Actividades planificadas por Proyecto.

% Parte 1 -> INICIALIZACIÓN

 ARP(i,:)=0;
 
  TR=0;
      
  k=0;
   
  SI=DataProject(:,4:6,i);
 
  ASASched=SI(ASched,:);                   % Act-> Conjunto de actividades  sucesoras inmediatas de las actividades del conjunto ASched j    
 
  Act=nonzeros(ASASched(1:numel(ASASched)))';
  
 TipoAct=(DataProject(Act,2,i))';     %TipoAct-> Tipos de las actividades del 

 WAct=Schedule(i,Act,4);                      % WAct -> Prioridades del conjunto de activiades Act

 OrdenAct=sortrows([Act;TipoAct;WAct]',3);   % OrdenAct-> Actividades Ordenadas segun WAct
 
 Act=OrdenAct(:,1)';

 ModeProc=DataProject(:,7:9,i);               %  ModeProc-> guarda la duración de las actividades según la especilización del recurso 

 aux=length(Act);                             % aux -> Cantidad de actividades de conjunto Act
 
 
 while k < aux
     
     TFS=Schedule(:,:,2);
     
     k=k+1;
     
     [AA,p]=find(DataProject(:,[4 5 6],i)==Act(k)); % Conjunto de actividades antesesoras de la actividad Act(k) -> antecesora de j 
          
      CAP=sum(APP(i,AA));                           % Calcula la cantidad de actividades antecesora ya planificadas de
      
      if CAP==length(AA) && APP(i,Act(k))==0        % si la actividad Act(k) está lista y nunca ha sido planificada
           
          ARP(i,Act(k))=Act(k);                     % Se inserta Act(k) en las actividades planificadas para se usada de indice
            
           APP(i,Act(k))=1;                         % Se archiva Act(k) como planificada
          
            TipoR=DataProject(Act(k),2,i);          % Definir el tipo  de recurso que requiere la actividad j del proyecto i 
           
            if TipoR==0 % Si la actividad j require un tipo de recurso=0 => última actividad del proyecto 
       
            Schedule(i,Act(k),1)=max(Schedule(i,AA',2));       % Tiempo de inicio de la actividad j=CA
               
            Schedule(i,Act(k),2)=Schedule(i,Act(k),1);         %  Tiempo de finalización de la actidad j=CA 
       
            else
                
                if and(k==1,isempty(RDO)==1)==1 || TR~=TipoR || isempty(RDO)==1
              
                TR=TipoR;
           
                [RDO]=RrHhAsig(DataProject,Schedule,RRHH,TipoR,AA,TFS,i);
               
                end
                                
                TFAmr=TFS(Schedule(:,:,3)==RDO(1));
                                                               
                if isempty(TFAmr)==1
                    
                    % Asignar Recursos a la actividad Act(k)
                    
                   Schedule(i,Act(k),3)=RDO(1);
                                                                        
                   Schedule(i,Act(k),1)=max(Schedule(i,AA,2));
                         
                   Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RDO(1),3));
                   
                   RDO=RDO(2:length(RDO));                            
                                   
                else
                     Schedule(i,Act(k),3)=RDO(1);
                    
                     Schedule(i,Act(k),1)=max(max(TFAmr),max(Schedule(i,AA',2))); % Calcula tiempo de inicio de la actividad j del proyecto i
            
                     Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RDO(1),3)); % 
                     
                     RDO=RDO(2:length(RDO));
                     
                end
                                            
            end
           
      end
     
 end

end

function [RDO]=RrHhAsig(DataProject,Schedule,RRHH,TipoR,AA,TFS,i)

CTypeR=RRHH(RRHH(:,2)==TipoR,1)'; % Conjunto de recursos  igual a TipoR de recurso que requiere la actividad que esta lista para procesarce

TIS=Schedule(:,:,1);

% Variables

cr=length(CTypeR); % Cantidad de recurso CTypeR

TIrW=0;            % Tiempo de inicio de cada vez que el recurso r ha sido asignado 

TFrW=0;            % Tiempo de fin de cada vez que el recurso r ha sido asignado 

FL=zeros(1,cr);    % Fecha de liberación de los recursos   
                                       
CoefA=FL;          % Coeficiente de asignación de cada recurso r

for r=CTypeR % Para todos los recursos de TipoR

    LastAA=max(TFS(i,AA));      %  Mayor tiempo de finalización de las actividades antecesoras de la actividad Act(k) 

    if isempty(nonzeros(Schedule(:,:,3)==r))==1
        
        FL(r)=0;
    else
                
        TIrW=TIS(Schedule(:,:,3)==r);   %  Tiempos de inicio de las actividades que se les ha asignado el recurso r
    
        TFrW=TFS(Schedule(:,:,3)==r);   %  Tiempos de finalización de las actividades que se les ha asignado el recurso r
              
        FL(r)=max(TFrW);                % Fecha de liberacion del recurso r 
    
    end

    switch RRHH(r,3)                               % Segun la especilizacion del recurso r disponible
            
            case 1                                 % si el recurso r tiene especilizacion 1 
               
                CE=0.1;                            % ponderacion de la especilizacion del recurso r 
                
                ATipoR=(DataProject(:,2,i)==TipoR);% Conjunto de las actividades de tipo (TipoR) del proyecto i (Ver la varainte de poner de todos los proyectos)
                
                TwT=sum(DataProject(ATipoR,7,i));  % Total de contenido de trabajo de las actividades de las actividades de tipo (TipoR) del proyecto i  
            
                SW=TFrW-TIrW;
        
                CW=(sum(SW)/TwT);               % Proporción del trabajo realizado por el recurso r 
        
        case 2
            
                CE=0.2;
                
                ATipoR=(DataProject(:,2,i)==TipoR); % Conjunto de las actiidades de tipo (TipoR) del proyecto i (Ver la varainte de poner de todos los proyectos)
                
                TwT=sum(DataProject(ATipoR,8,i));   % Total de contenido de trabajo de las actividades de las actividades de tipo (TipoR) del proyecto i  
                
                SW=TFrW-TIrW;
        
                CW=(sum(SW)/TwT);               % Proporción del trabajo realizado por el recurso r 
        
        case 3
            
                CE=0.3;
                
                 ATipoR=(DataProject(:,2,i)==TipoR);% Conjunto de las actiidades de tipo (TipoR) del proyecto i (Ver la varainte de poner de todos los proyectos)
                
                TwT=sum(DataProject(ATipoR,9,i));   % Total de contenido de trabajo de las actividades de las actividades de tipo (TipoR) del proyecto i  
                
                 SW=TFrW-TIrW;
        
                CW=(sum(SW)/TwT);               % Proporción del trabajo realizado por el recurso r 
     end

         if FL(r)==0 && max(max(TFS))==0
                         
            CL=0;
       
         else
             
            CL=(1-(max(max(TFS))-FL(r))/max(max(TFS)));
            
         end

         if FL(r) > LastAA      % Si el recurso r esta disponible
             
             ND=0.5;            % Ponderación del recurso r si no está disponible
             
         else 
               
             ND=0;               % Ponderación del recurso r si está disponible
       
         end
                              
        CoefA(CTypeR==r)=(CE+CL+CW+ND); % CE->ponderacíón de Tipo (especilización) Fecha de Liberación
        
end

    aux1=sortrows([CoefA;CTypeR]',1); % Recursos ordenados segun Coeficiente de asignacion 

    RDO=aux1(:,2)';
    
end

function [Schedule]=NeighborhoodR(DataProject,Schedule,NP,CA)
% NeighborhoodR function  
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
% Facultad III. Universidad de las Ciencias Informáticas
% Descripción: Realiza una intercambio de recursos del mismo tipo entre dos
% actividades cualquiera.
% ENTRADAS
%
% 
% IMPLEMENTACIÓN 

aux=1:NP;

SelectAct1=unidrnd(CA);

SelectProject1=unidrnd(NP);

if SelectAct1==1 || SelectAct1==CA  
    
return

else
            
    TipoA=DataProject(SelectAct1,2,SelectProject1); % Tipo de actividad = tipo de recurso que requiere en la primera selección
    
    % [CATipoA,u]=find(DataProject(:,2,SelectProject1)==TipoA);   % Conjunto de las actividades TipoA en  la primera selección
    
    aux(SelectProject1)=0;
    
    aux2=nonzeros(aux)';
    
    SelectProject2=aux2(unidrnd(length(aux2)));
    
   [CATipoA2,u]=find(DataProject(:,2,SelectProject2)==TipoA); % conjunto de actividades de TipoA en la Segunda selección
   
    SelectAct2=CATipoA2(unidrnd(length(CATipoA2)));
    
    RAct2=Schedule(SelectProject2,SelectAct2,3);
    
    if RAct2~=Schedule(SelectProject1,SelectAct1,3)
        
        aux3=Schedule(SelectProject1,SelectAct1,3);
        
        Schedule(SelectProject1,SelectAct1,3)=RAct2;
        
        Schedule(SelectProject2,SelectAct2,3)=aux3;
        
    else
        
        return
        
    end 
    
    
end
end

function [Schedule,FO] = SSGSSA(Schedule,DataProject,RRHH,TA,CP,NP,CA)
% SSGSnew Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
%  Facultad III. Universidad de las Ciencias Informáticas
%  Descripción: Heuristica para la generación de soluciones iniciales factibles basada
% en: % Serial Scheme Generation Sechdeule SSGS  
% Diseñada para  n proyectos n=2,3,... con la misma cantidad de actividades
% con tiempos de arribo Start=[T1,T2,T3,T4,...,Tn]
% los proyectos son tomados con dos actividades ficticias
% ENTRADA
%  DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Información sobre los recursos 
%  TA -> Tiempos de arribo de los proyectos  TA=[TA1,TA2,...,TAn]
% SALIDAS
%  Schedule -> Solución  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  TI -> Tiempos de Inicio
%  TF -> Tiempos de Finalización 
%  RA -> Recurso Asignado

% DEFINICIÓN DE VARIABLES

Schedule(:,:,1)=0;

Schedule(:,:,2)=0;

PS=zeros(1,NP); % Variable que guarda los proyectos que han comenzado PS(i)={0,1} 

PF=PS;      % Variable que guarda los proyectos que han finalizado PF(i)={0,1} 

APP=zeros(NP,CA);  % Variable que guarda las actividades que se han planificado APP(i,j)={0,1} 

% INICIALIZACIÓN

PS(1)=1; %  Indica que el proyecto 1 comienza procesandose

APP(1,1)=1; %  Se planifica la primera actividad del primer proyecto

ARP=APP; % Conjuntos de actividades que pueden ser planificadas se inicializa con  la actividad 1 del proyecto 1  
 
Stop=0; %  Se inicializa en cero la condición de parada

ind=2:NP; % indice de los proyectos 

ind2=1:NP;

v=1; % Contador de los proyectos que inician

while Stop < NP % Mientras la cantidad de proyectos finalizado sea menor que la actidad de proyectos 
    
if v~=NP 

    MK=max(max(Schedule(:,2:end,2)));
            
for k=nonzeros(ind)' % Para cada uno de los Proyectos
    
    if TA(k)<=MK && TA(k)~=0 % Si el mayor si el Makespan parcial MK es mayor o igual a la fecha de arribo deque la fe
          
       PS(k)=1;
       
       Schedule(k,1,1)=TA(k); % se inician los proyectos que aun no han comenzado
       
       Schedule(k,1,2)=TA(k);
              
       v=v+1;
       
       ind(k-1)=0;
       
       ARP(k,1)=1;
       
       APP(k,1)=1;
                    
     end
       
end
    
end
       
    PR=find(PS==1);% Proyectos que han comenzado 
    
    CPS=length(PR); % Cantidad de proyectos iniciados
   
    if CPS==1 % Si la cantidad de proyectos iniciados es 1 
         
         OrdenP=PR; % Orden es PR
         
    else % en caso contrario
        
         OrdenP=zeros(1,length(PR));
          
         aux=rand(1,CPS); % generar CPS numeros aleatorios
         
         ORD=sort(aux);   %  Ordenar 
         
         for l=1:length(PR)
                     
              OrdenP(l)=PR(aux==ORD(l)); % Establecer orden
                            
         end
         
    end

    for i=OrdenP % Para proyectos Iniciado segun OrdenP
                              
     [Schedule,APP,ARP] = SchedulingFeasibleSA(Schedule,DataProject,RRHH,i,nonzeros(ARP(i,:))',ARP,APP);                                                         
                                              
    end
        
  for r=nonzeros(ind2)'
    
    if Schedule(r,CA,2)~=0 && sum(APP(r,:))==CA
        
        PF(r)=1;
        
        PS(r)=0;   
        
        ind2(r)=0;
        
        ARP(r,:)=0;
             
    end
    
  end
  Stop=sum(PF);
end  
FO=mean(((Schedule(:,CA,2)'-(CP+TA))));  %mean((Schedule(:,CA,2)'./CP)-1)); %mean(Schedule(:,CA,2));
end

function [Schedule,APP,ARP] = SchedulingFeasibleSA(Schedule,DataProject,RRHH,i,j,ARP,APP)
%  SchedulingFeasible Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias Básicas.
%  Facultad III. Universidad de las Ciencias Informáticas
%  Planifica las activitades de ActReadyPlan sin violar las restricciones de
% Precedencia y de Recursos.
% ERNTRADAS
% Schedule -> solución 
% DataProyect -> Datos de los Proyectos 
%ActReadyPlan-> Conjunto de actividades que están lista para ser planificadas 
% RRHH-> Datos sobre los recursos humanos del proyecto
% APP(i) -> Conjunto de Actividades planificadas por Proyecto.
% DataProject -> Informaci�n sobre los proyectos.
% SALIDAS
% Schedule-> solución 
% APP(i) -> Actividades planificadas por Proyecto.

% INICIALIZACIÓN

 ARP(i,:)=0;
 
 AR=nonzeros(DataProject(nonzeros(j),4:6,i))'; % Conjunto de actividades  sucesoras inmediatas de la actividad j    
 
 ModeProc=DataProject(:,7:9,i);

 aux=length(AR);
 
 AxW=[Schedule(i,AR,4);AR];
 
   OrdenA=sortrows(AxW',1);
 
   Act=OrdenA(:,2)';
 
 for k=1:aux
     
      [AA,p]=find(DataProject(:,[4 5 6],i)==Act(k));
          
      CAP=sum(APP(i,AA'));
      
      if CAP==length(AA) && APP(i,Act(k))==0 % si está lista y nunca ha sido planificada
          
          ARP(i,Act(k))=Act(k);
            
           APP(i,Act(k))=1;
          
          TipoR=DataProject(Act(k),2,i); % Definir el tipo  de recurso que requiere la actividad j del proyecto i
           
          if TipoR==0 % Si la actividad j require un tipo de recurso=0 => última actividad del proyecto 
       
            Schedule(i,Act(k),1)=max(Schedule(i,AA',2)); % Tiempo de inicio de la actividad j=CA
               
            Schedule(i,Act(k),2)=Schedule(i,Act(k),1);         %  Tiempo de finalización de la actidad j=CA 
                      
          else
               
             RAsig=Schedule(i,Act(k),3); % Asignación del Recurso
            
             TFS=Schedule(:,:,2);
                                       
             TFAmr=TFS(Schedule(:,:,3)==RAsig); % Buscar si el recurso RAsig ha sido asignado anteriormente
             
             if isempty(TFAmr)==1 && RAsig~=0 % Si el recurso asignado a la actividad j nunca ha sido asignado
                 
               Schedule(i,Act(k),1)=max(Schedule(i,AA',2)); % Calcula tiempo de inicio de la actividad j del proyecto i
               
                             
               Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalización de la actividad j del proyecto i
             
             else
                if RAsig~=0
                                 
               Schedule(i,Act(k),1)=max(max(TFAmr),max(Schedule(i,AA',2))); % Calcula tiempo de inicio de la actividad j del proyecto i
             
               Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalización de la actividad j del proyecto i
                end
             end
          end          
      end
 end
end
