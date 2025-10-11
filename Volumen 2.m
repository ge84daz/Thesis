%�bungsbeispiel Volumenelemente
%by Martin Buchschmid und Siegfried Seipelt

clear all;

femesh('reset');

%Deklaration der ben�tigten Startknoten
FEnode=[1  0 0 0  0 0 0 ;
        2  0 0 0  0 0.2 0 ;
        3  0 0 0  0 1 0];%nicht verwendete Knoten stellen kein Problem dar
        
%Vorgabe der zun�chst leeren Felder f�r die FE-Elementierung 
FEelt=[];
FEel0=[];

%Erzeugen eines Balkenelementes zwischen den Knoten 1 und 2
femesh('objectbeamline 1 2');

%Extrudieren des Grundelementes jeweils 10-fach in die Richtungen x=1 und
%dann z=1 (Die zuerst erzeugte Fl�che wird so zu einem Volumenelement)
femesh('extrude 30 0.2 0 0');
femesh('extrude 30 0 0 0.2');

%Anf�gen der zuletzt verwendeten Elemente in das Arbeitsfeld
femesh('addsel');

%Verschieben der Volumenelementscheibe um y=1 (mit transsel 10 0 1 0 w�rde
%die selbe Aktion 10 mal ausgef�hrt werden
femesh('transel 0 0.2 0');

% Hinzuf�gen der Geometriedaten zum endg�ltigen Feld
femesh('addsel');

%�bergabe der Geometriedaten des FE-Netzes in das f�r die Berechnung
%verwendete Feldt [model]
model=femesh('model');

%Eingabe der Materialkennwerte
model.pl=m_elastic('dbval 1 Steel');

% Bei 3D-Elementen erfolgt hier eine vorgabe der Integrationsregel, welche
% vom gew�hlten Elementtyp abh�ngt
model.il=[1 fe_mat('p_solid','SI',1) 0 3 0 2];  

load=struct('DOF',336.02,'def',50);

% Eingabe der Randbedingungen 
model = fe_case(model,'fixdof','linker Rand eingespannt','x==0',...
                      'DofLoad','Punktlast',load);
def1 = fe_eig(model);
cf1=feplot(model,def1);
% Wahl der Time-Step-Methode und Vorgabe der zugeh�rigen Parameter
  TimeOpt=struct('Method','Newmark','Opt',[.25 .5 0 .00001 300 5],'NeedUVA',[1 0 0]);
  %                        a               b   c  d e    f   g   h        i
  % a: Mehtode
  % b: alfa-Faktor
  % c: beta-Faktor
  % d: Zeit der Lastaufbringung in Sekunden
  % e: Zeitschrittbreite
  % f: Anzahl der Zeitschritte
  % g: Einwirkungsdauer der transienten Last
  % h: Vorgabe der Ausgabewerte U=Verformung, V=Geschwindigkeit, A=Beschl.
  % i: 1 bedeutet Ausgabe, 0 bedeutet keine Ausgabe 
  
  % Assemblierung des Models mit der gew�hlten Methode und 
  % Speicherung der Daten in der .def Matrix
  def=fe_time(TimeOpt,model);
    
  % Plotten des Models und Animation
  cf=feplot;cf.model=model;cf.def=def;
  fecom(';view1;animtime;ch1');
  % ch20 bedeutet, dass die Animation beim 20 Timestep beginnt
