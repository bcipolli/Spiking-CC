#/bin/bash
rdir=izzy/time
echo "put spnet.cpp $rdir/spnet.cpp" | sftp dino.ucsd.edu
ssh dino.ucsd.edu "g++ -o $rdir/spnet $rdir/spnet.cpp"
echo "get $rdir/spnet" | sftp dino.ucsd.edu
