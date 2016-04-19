
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <streambuf>

mysint rnk(mysint randNk){
    mysint random_integer=rand()%randNk;
    return random_integer;
}
myint getNextInteger(char * memblock, std::streampos size, std::streampos *pos){
    
    myint current=0;
    myint mult=1;
    
    while ((*pos<size) &&(memblock[*pos]!='-')&&(( memblock[*pos]<'0')||(memblock[*pos]>'9'))){
        *pos= (*pos)+((std::streampos)1);
    }
    if(memblock[*pos]=='-'){
        mult=-1;
        *pos= (*pos)+((std::streampos)1);
    }
    while ((*pos<size) &&(( memblock[*pos]>='0')&&(memblock[*pos]<='9'))){
        current= 10 * current;
        current+= (int)(memblock[*pos]-'0');
        *pos=*pos+((std::streampos)1);
    }
    current=current*mult;
    return current;
}


myint readSequenceFromFile(std::string filename, mysint** sequence1, myint *length){
    srand((unsigned)time(0));
    std::streampos size; std::streampos *position;
    position=new std::streampos;
    *position=0;
    char * memblock;
    std::ifstream ifile(filename,std::ios::in|std::ios::binary|std::ios::ate);
    if (ifile.is_open())
    {
        size = ifile.tellg();
        memblock = new char [size];
        ifile.seekg (0, std::ios::beg);
        ifile.read (memblock, size);
        ifile.close();
        *length=getNextInteger(memblock,size,position);
        *sequence1=new mysint[*length];
        myint i=0;mysint cNum=0;
        while((i<*length)&&(cNum!=-9)){
            cNum=getNextInteger(memblock,size,position);
            if(cNum!=-9){
                (*sequence1)[i]=cNum;
                i++;
            }
        }
        
        while(i<(*length)){
            (*sequence1)[i]=rnk(20);
            i++;
            
        }
        
        delete[] memblock;
        
    }
    return 1;
}


myint printToFile(std::string filename, myint *s, myint l){
    
    std::ofstream mfile;
    mfile.open(filename);
    mfile << l<<std::endl;
    for(myint i=0;i<l;i++){
        mfile<<s[i]<<" ";
    }
    mfile<<std::endl;
    mfile.close();
    return 1;
}
