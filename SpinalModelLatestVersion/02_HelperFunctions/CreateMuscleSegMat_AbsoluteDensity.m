%% Build Muscle Matrix

Muscles = {'Infrahyoid','Sternomastoid','Axial','Trapeze','Phrenic','Rhomboid','Deltoid','Supraspinatus','Biceps','Serratus','Pectoral','Triceps','Latissmus Dorsi','Forearm Flex','Forearm Ext','Intercostal','Thoracic Abdominal',...
'Quadratus Lumburum','Rectus Femoris','Vastus Lateralis','Vastus Medialis','Iliopsoas','Adductors','Gracilis','Semitendinosus','Biceps Femoris','Gluteus','Tibialis Anterior','Soleus','Gastrocnemius','Extensor Digitorum','Plantar Flexor','Tail'};
FlexExt = categorical({'','','','','','','','','Flex','','','Ext','','Flex','Ext','','','Flex','Bi','Ext','Ext','Flex','Flex','Bi','Bi','Flex','Ext','Flex','Ext','Ext','Ext','Ext',''});
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

T2.T13('Axial') = 10;
T2.T13('Thoracic Abdominal') = 10;
T2.T13('Vastus Lateralis') = 5;

T2.L1('Quadratus Lumburum') = 10;
T2.L1('Axial') = 10;
T2.L1('Vastus Lateralis') = 10;
T2.L1('Vastus Medialis') = 5;

T2.L2('Rectus Femoris') = 10;
T2.L2('Vastus Lateralis') = 15;
T2.L2('Vastus Medialis') = 15;
T2.L2('Adductors') = 7;
T2.L2('Iliopsoas') = 40;
T2.L2('Axial') = 10;

T2.L3('Rectus Femoris') = 35;
T2.L3('Vastus Medialis') = 35;
T2.L3('Vastus Lateralis') = 35;
T2.L3('Iliopsoas') = 15;
T2.L3('Gastrocnemius') = 5;
T2.L3('Adductors') = 36;
T2.L3('Gracilis') = 29;
T2.L3('Gluteus') = 11;
T2.L3('Axial') = 10;
T2.L3('Semitendinosus') = 20;
T2.L3('Tibialis Anterior') = 7;


T2.L4('Rectus Femoris') = 5;
T2.L4('Vastus Medialis') = 5;
T2.L4('Semitendinosus') = 28;
T2.L4('Biceps Femoris') = 15;
T2.L4('Axial') = 10;
T2.L4('Gluteus') = 21;
T2.L4('Tibialis Anterior') = 50;
T2.L4('Soleus') = 5;
T2.L4('Gastrocnemius') = 15;
T2.L4('Adductors') = 18;
T2.L4('Gracilis') = 16;
T2.L4('Iliopsoas') = 5;
T2.L4('Plantar Flexor') = 5;

T2.L5('Soleus') = 40;
T2.L5('Gastrocnemius') = 35;
T2.L5('Semitendinosus') = 24;
T2.L5('Biceps Femoris') = 33;
T2.L5('Gluteus') = 34;
T2.L5('Adductors') = 10;
T2.L5('Axial') = 10;
T2.L5('Extensor Digitorum') = 4;
T2.L5('Plantar Flexor') = 20;
T2.L5('Gracilis') = 8;
T2.L5('Tibialis Anterior') = 20;

T2.L6('Soleus') = 30;
T2.L6('Gastrocnemius') = 15;
T2.L6('Gluteus') = 10;
T2.L6('Semitendinosus') = 17;
T2.L6('Biceps Femoris') = 14;
T2.L6('Extensor Digitorum') = 45;
T2.L6('Plantar Flexor') = 10;

T2.S1('Tail') = 10;
T2.S2('Tail') = 10;
T2.S3('Tail') = 10;
T2.S4('Tail') = 10;

MuscleSeg = T2;
%MuscleSeg.Variables = normalize(MuscleSeg.Variables,2,'norm');%./sum(MuscleSeg.Variables,1);
save(".\01_SpinalGeometryFiles\matMuscleSeg.mat","MuscleSeg","FlexExt");

