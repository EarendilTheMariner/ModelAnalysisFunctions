%% Build Muscle Matrix

Muscles = {'Infrahyoid','Sternomastoid','Axial','Trapeze','Phrenic','Rhomboid','Deltoid','Supraspinatus','Biceps','Serratus','Pectoral','Triceps','Latissmus Dorsi','Forearm Flex','Forearm Ext','Intercostal','Thoracic Abdominal',...
'Quadratus Lumburum','Rectus Femoris','Vastus Lateralis','Vastus Medialis','Iliopsoas','Adductors','Gracilis','Semitendinosus','Biceps Femoris','Gluteus','Tibialis Anterior','Soleus','Gastrocnemius','Extensor Digitorum','Plantar Flexor','Tail'};
FlexExt = categorical({'','','','','','','','','Flex','','','Ext','','Flex','Ext','','','Flex','Bi','Ext','Ext','Flex','Flex','Bi','Bi','Flex','Ext','Flex','Ext','Ext','Ext','Flex',''});
load("01_SpinalGeometryFiles\matSegmentDimensions.mat")

T2 = table('Size',[length(Muscles),length(Segments)],'VariableTypes',repmat("double",[length(Segments),1]),'VariableNames',Segments,'RowNames',Muscles);

T2.C1('Infrahyoid') = 0.3;
T2.C1('Axial') = 1;
T2.C1('Sternomastoid') = 1;

T2.C2('Infrahyoid') = 1;
T2.C2('Trapeze') = 0.1;
T2.C2('Axial') = 1;

T2.C3('Axial') = 1;
T2.C3('Phrenic') = 0.3;
T2.C3('Trapeze') = 1;
T2.C3('Infrahyoid') = 0.25;

T2.C4('Deltoid') = 0.3;
T2.C4('Trapeze') = 0.1;
T2.C4('Phrenic') = 1;
T2.C4('Rhomboid') = 0.5;
T2.C4('Axial') = 1;

T2.C5('Deltoid') = 1;
T2.C5('Biceps') = 0.5;
T2.C5('Supraspinatus') = 0.5;
T2.C5('Phrenic') = 0.1;
T2.C5('Axial') = 1;
T2.C5('Rhomboid') = 1;

T2.C6('Deltoid') = 0.3;
T2.C6('Supraspinatus') = 0.5;
T2.C6('Biceps') = 0.5;
T2.C6('Serratus') = 0.50;
T2.C6('Axial') = 1;

T2.C7('Axial') = 1;
T2.C7('Serratus') = 0.5;
T2.C7('Pectoral') = 0.5;
T2.C7('Triceps') = 0.5;
T2.C7('Latissmus Dorsi') = 0.7;
T2.C7('Forearm Flex') = 0.8;
T2.C7('Forearm Ext') = 0.2;

T2.C8('Forearm Ext') = 0.2;
T2.C8('Forearm Flex') = 0.8;
T2.C8('Latissmus Dorsi') = 0.3;
T2.C8('Triceps') = 1;
T2.C8('Pectoral') = 1;
T2.C8('Axial') = 1;

T2.T1('Axial') = 1;
T2.T1('Triceps') = 0.3;
T2.T1('Pectoral') = 0.3;
T2.T2('Axial') = 1;
T2.T2('Intercostal') = 1;
T2.T3 = T2.T2;
T2.T4 = T2.T2;
T2.T5 = T2.T2;
T2.T6 = T2.T2;
T2.T7 = T2.T2;
T2.T8 = T2.T2;
T2.T9('Axial') = 1;
T2.T9('Thoracic Abdominal') = 1;
T2.T10 = T2.T9;
T2.T11 = T2.T9;
T2.T12('Thoracic Abdominal') = 1;
T2.T12('Axial') = 1;

T2.T13('Axial') = 1;
T2.T13('Thoracic Abdominal') = 0.3;

T2.L1('Quadratus Lumburum') = 1;
T2.L1('Axial') = 1;
T2.L1('Rectus Femoris') = 0.08;
T2.L1('Vastus Lateralis') = 0.08;
T2.L1('Vastus Medialis') = 0.08;

T2.L2('Rectus Femoris') = 0.085;
T2.L2('Vastus Lateralis') = 0.085;
T2.L2('Vastus Medialis') = 0.085;
T2.L2('Adductors') = 0.093;
T2.L2('Iliopsoas') = 0.17;
T2.L2('Axial') = 1;

T2.L3('Rectus Femoris') = 0.68;
T2.L3('Vastus Medialis') = 0.68;
T2.L3('Vastus Lateralis') =0.68;
T2.L3('Iliopsoas') = 0.73;
T2.L3('Gastrocnemius') = 0.2;
T2.L3('Adductors') = 0.5;
T2.L3('Gracilis') = 0.54;
T2.L3('Gluteus') = 0.15;
T2.L3('Axial') = 1;
T2.L3('Semitendinosus') = 0.18;

T2.L4('Rectus Femoris') = 0.15;
T2.L4('Vastus Medialis') = 0.15;
T2.L4('Vastus Lateralis') =0.15;
T2.L4('Semitendinosus') = 0.33;
T2.L4('Biceps Femoris') = 0.14;
T2.L4('Axial') =1;
T2.L4('Gluteus') = 0.3;
T2.L4('Tibialis Anterior') = 0.53;
T2.L4('Soleus') = 0.08;
T2.L4('Gastrocnemius') = 0.3;
T2.L4('Adductors') = 0.25;
T2.L4('Gracilis') = 0.3;
T2.L4('Iliopsoas') = 0.09;
T2.L4('Plantar Flexor') = 0.07;

T2.L5('Soleus') = 0.55;
T2.L5('Gastrocnemius') = 0.5;
T2.L5('Semitendinosus') = 0.28;
T2.L5('Biceps Femoris') = 0.6;
T2.L5('Gluteus') = 0.48;
T2.L5('Adductors') = 0.15;
T2.L5('Axial') = 0.5;
T2.L5('Extensor Digitorum') = 0.05;
T2.L5('Plantar Flexor') = 0.6;
T2.L5('Gracilis') = 0.16;
T2.L5('Tibialis Anterior') = 0.47;

T2.L6('Soleus') = 0.37;
T2.L6('Gluteus') = 0.08;
T2.L6('Semitendinosus') = 0.2;
T2.L6('Biceps Femoris') = 0.26;
T2.L6('Extensor Digitorum') = 0.95;
T2.L6('Plantar Flexor') = 0.3;

T2.S1('Tail') = 1;
T2.S2('Tail') = 1;
T2.S3('Tail') = 1;
T2.S4('Tail') = 1;

MuscleSeg = T2;
MuscleSeg.Variables = MuscleSeg.Variables./sum(MuscleSeg.Variables,2);
save(".\01_SpinalGeometryFiles\matMuscleSeg.mat","MuscleSeg","FlexExt");

