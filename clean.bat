rd .vs /S /Q

cd "1. CPU Path Tracer Vanilla"
rd x64\debug /S /Q
rd x64\release /S /Q
rd x64 /S /Q
del buildlog.txt
cd ..

cd "2. CPU Path Tracer Enhanced"
rd x64\debug /S /Q
rd x64\release /S /Q
rd x64 /S /Q
del buildlog.txt
cd ..

cd "3. GPU Path Tracer"
rd x64\debug /S /Q
rd x64\release /S /Q
rd x64 /S /Q
del buildlog.txt
cd ..