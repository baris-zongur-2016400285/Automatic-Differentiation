#include <iostream>
#include <math.h>
#include <vector>
#include <sstream>
#include <fstream>
using namespace std;
class Node{
public:
    //I didnt used template
    //Each node keeps a function by its name , required other vars or funcs fronts and backs;
    string name;
    double myVal;
    double diffVal;
    string fnc="nof";
    int indegr;
    vector<double> vars;
    vector<string> myVals;
    vector<Node*> fronts;
    vector<Node*> backs;
    bool iscalc=false;
    vector<string> incvals;
    Node(string name){
        this->name = name;
    };
    void calculate(){
        //this part calculates the value of functions
        if(fnc.compare("nof")!=false){
            if(fnc.compare("mult")==0){
                myVal = vars[0]*vars[1];
            }

            if(fnc.compare("add")==0){
                myVal = vars[0]+vars[1];
            }

            if(fnc.compare("subs")==0){
                myVal = vars[0]-vars[1];
            }

            if(fnc.compare("divide")==0){
                myVal = vars[0]/vars[1];
            }

            if(fnc.compare("sin")==0){
                myVal=sin(vars[0]);
            }
            if(fnc.compare("cos")==0){
                myVal=cos(vars[0]);
            }
            if(fnc.compare("tan")==0){
                myVal=tan(vars[0]);
            }
            if(fnc.compare("acos")==0){
                myVal=acos(vars[0]);
            }
            if(fnc.compare("asin")==0){
                myVal=asin(vars[0]);
            }
            if(fnc.compare("atan")==0){
                myVal=atan(vars[0]);
            }
            if(fnc.compare("exp")==0) {
                myVal=exp(vars[0]);
            }
            if(fnc.compare("log")==0){
                myVal=log(vars[0]);
            }
            if(fnc.compare("log10")==0){
                myVal=log10(vars[0]);
            }
            if(fnc.compare("pow")==0){
                myVal=pow(vars[0],vars[1]);
            }
            if(fnc.compare("sqrt")==0){
                myVal=sqrt(vars[0]);
            }
        }
    }
    double calculateDiff(string myinp){
        //this part calculates diff wrt which var i want
        if(fnc.compare("nof")==0){
            if(name.compare(myinp)==0){
                return 1.0;
            }
            else{
                return 0;
            }
        }
            if (fnc.compare("nof") != 0) {
                if (fnc.compare("mult") == 0) {
                   return (backs[0]->calculateDiff(myinp)*vars[1])+(backs[1]->calculateDiff(myinp)*vars[0]);
                }

                if (fnc.compare("add") == 0) {
                  return (backs[0]->calculateDiff(myinp))+(backs[1]->calculateDiff(myinp));
                }

                if (fnc.compare("subs") == 0) {
                    return (backs[0]->calculateDiff(myinp))-(backs[1]->calculateDiff(myinp));
                }

                if (fnc.compare("divide") == 0) {
                    return ((backs[0]->calculateDiff(myinp)*vars[1])-(backs[1]->calculateDiff(myinp)*vars[0]))/((vars[1]*vars[1]));
                }

                if (fnc.compare("sin") == 0) {
                   return cos(vars[0])*backs[0]->calculateDiff(myinp);
                }
                if (fnc.compare("cos") == 0) {
                    return sin(vars[0])*backs[0]->calculateDiff(myinp)*(-1);
                }
                if (fnc.compare("tan") == 0) {
                   return backs[0]->calculateDiff(myinp)*(1+pow(tan(vars[0]),2));
                }
                if (fnc.compare("acos") == 0) {
                    return ((-1)*backs[0]->calculateDiff(myinp))/(sqrt(1-pow(vars[0],2)));
                }
                if (fnc.compare("asin") == 0) {
                    return backs[0]->calculateDiff(myinp)/(sqrt(1-pow(vars[0],2)));
                }
                if (fnc.compare("atan") == 0) {
                   return backs[0]->calculateDiff(myinp)/(1+pow(vars[0],2));
                }
                if (fnc.compare("exp") == 0) {
                    return backs[0]->calculateDiff(myinp)*exp(vars[0]);
                }
                if (fnc.compare("log") == 0) {
                    return backs[0]->calculateDiff(myinp)/vars[0];
                }
                if (fnc.compare("log10") == 0) {
                    return backs[0]->calculateDiff(myinp)/(vars[0]*log(10));
                }
                if (fnc.compare("pow") == 0) {
                   return pow(vars[0],vars[1])*(backs[1]->calculateDiff(myinp)*log(vars[0])+(backs[0]->calculateDiff(myinp)/vars[0])*vars[1]);
                }
                if (fnc.compare("sqrt") == 0) {
                    return (0.5)*(1.0/sqrt(vars[0]))*backs[0]->calculateDiff(myinp);
                }
            }
        }
};
class Graph{
public:
    vector<Node*> nodes;
    string myfunc;
    vector<string> variables;
    bool topoSort(){
        //this is a classic topo sort algorith turns true if no cycly
        // turns false if there is a cycle
        vector<Node*> topoSort;
        int cnt =0;
        for (int i=0;i<nodes.size();i++){
            nodes[i]->indegr=nodes[i]->backs.size();
            if(nodes[i]->indegr==0){
                topoSort.push_back(nodes[i]);
                cnt++;
            }
        }
        while(topoSort.empty()==false){
            for(int k=0;k<topoSort.size();k++){
                if(topoSort[k]->indegr==0){
                    for(int l=0;l<topoSort[k]->fronts.size();l++){
                        topoSort[k]->fronts[l]->indegr--;
                        if(topoSort[k]->fronts[l]->indegr==0){
                            topoSort.push_back(topoSort[k]->fronts[l]);
                            cnt++;
                        }

                    }
                topoSort.erase(topoSort.begin()+k);
                    break;
                }

            }

        }
        if(cnt==nodes.size()){
            return true;
        }
        else{
            return false;
        }
    }
};


int main(int argc, char* argv[]) {
    vector<string> inputS;
    vector<string> inputS2;
    string dervs="";
    string result="";
    int iscyc=0;
    if(argc != 5){ //this part reads the input and checks if it is ok
        cout<<"terminal error"<<endl;
    }
    else {
        ifstream infile(argv[1]); //classic read line by line
        if (!infile.is_open()) {
            cout << "couldn't open file" << endl;
            return 0;
        } else {
            string s;
            while (getline(infile, s)) {
                inputS.push_back(s);
            }
        }
        infile.close();
    }
    if(argc != 5){ //this part reads the input and checks if it is ok
        cout<<"terminal error"<<endl;
    }
    else {
        ifstream infile(argv[2]); //classic read line by line
        if (!infile.is_open()) {
            cout << "couldn't open file" << endl;
            return 0;
        } else {
            string s;
            while (getline(infile, s)) {
                inputS2.push_back(s);
            }
        }
        infile.close();
    }
    vector<string> variables;
    string nowuse="";
    for(int i=0;i<inputS2[0].length();i++){
        if(inputS2[0][i]!=' '){
            nowuse = nowuse +inputS2[0][i];
        }
        else{
            variables.push_back(nowuse);
            nowuse="";
        }
    }
    for(int i=1;i<inputS2.size();i++){
        //for each different inputs creates graph and calculets required values
        Graph myGraph;
        vector<double> myvars;
        string nowuse="";
        for(int j=0;i<inputS2[i].length();j++){
            if(myvars.size()==variables.size()){
                break;
            }
            if(inputS2[i][j]!=' '){
                nowuse = nowuse +inputS2[i][j];
            }
            else{
                myvars.push_back(stod(nowuse));
                nowuse="";
            }
        }
        for(int i=0;i<variables.size();i++){
            Node* myNode = new Node(variables[i]);
            myGraph.nodes.push_back(myNode);
            myGraph.nodes[i]->myVal=myvars[i];
        }
        nowuse="";
        vector<string> nowvars;
        for(int i=0;i<inputS[variables.size()].length();i++){
            if(inputS[variables.size()][i]!=' '){
                nowuse = nowuse +inputS[variables.size()][i];
            }
            else{
                nowvars.push_back(nowuse);
                nowuse="";
            }
        }
        nowvars.push_back(nowuse);
        myGraph.myfunc=nowvars[1];
        nowvars.clear();
        nowuse="";
        for(int i = variables.size()+1;i<inputS.size();i++){
            nowvars.clear();
            nowuse="";
            for(int j=0;j<inputS[i].size();j++){
                if(inputS[i][j]!=' '){
                    nowuse = nowuse +inputS[i][j];
                }
                else{
                    nowvars.push_back(nowuse);
                    nowuse="";
                }
                if(inputS[i][j]=='\\'){
                    break;
                }
            }
            nowvars.push_back(nowuse);
            nowuse="";
            if(nowvars.size()==4){
                Node* myNode = new Node(nowvars[0]);
                myGraph.nodes.push_back(myNode);
                for(int k=0;k<myGraph.nodes.size();k++){
                    if(nowvars[3].compare(myGraph.nodes[k]->name)==0){
                        myGraph.nodes[myGraph.nodes.size()-1]->fnc=nowvars[2];
                        myGraph.nodes[myGraph.nodes.size()-1]->vars.push_back(myGraph.nodes[k]->myVal);
                        Node* perm =myGraph.nodes[k];
                        myGraph.nodes[myGraph.nodes.size()-1]->backs.push_back(perm);
                        myGraph.nodes[myGraph.nodes.size()-1]->myVals.push_back(myGraph.nodes[k]->name);
                        Node* perm2 =myGraph.nodes[myGraph.nodes.size()-1];
                        myGraph.nodes[k]->fronts.push_back(perm2);
                        myGraph.nodes[myGraph.nodes.size()-1]->calculate();
                        break;
                    }
                }
            }
            if(nowvars.size()==5){
                Node* myNode = new Node(nowvars[0]);
                myGraph.nodes.push_back(myNode);
                myGraph.nodes[myGraph.nodes.size()-1]->fnc=nowvars[2];
                for(int k=0;k<myGraph.nodes.size();k++){
                    if(nowvars[3].compare(myGraph.nodes[k]->name)==0){
                        myGraph.nodes[myGraph.nodes.size()-1]->vars.push_back(myGraph.nodes[k]->myVal);
                        Node* perm =myGraph.nodes[k];
                        myGraph.nodes[myGraph.nodes.size()-1]->backs.push_back(perm);
                        myGraph.nodes[myGraph.nodes.size()-1]->myVals.push_back(myGraph.nodes[k]->name);
                        Node* perm2 = myGraph.nodes[myGraph.nodes.size()-1];
                        myGraph.nodes[k]->fronts.push_back(perm2);
                        break;
                    }
                }
                for(int k=0;k<myGraph.nodes.size();k++){
                    if(nowvars[4].compare(myGraph.nodes[k]->name)==0){
                        myGraph.nodes[myGraph.nodes.size()-1]->vars.push_back(myGraph.nodes[k]->myVal);
                        Node* perm =myGraph.nodes[k];
                        myGraph.nodes[myGraph.nodes.size()-1]->backs.push_back(perm);
                        myGraph.nodes[myGraph.nodes.size()-1]->myVals.push_back(myGraph.nodes[k]->name);
                        Node* perm2 =myGraph.nodes[myGraph.nodes.size()-1];
                        myGraph.nodes[k]->fronts.push_back(perm2);
                        break;
                    }

                }
                myGraph.nodes[myGraph.nodes.size()-1]->calculate();
            }
        }
        //check if cyclic
        if(iscyc==0){
            if(myGraph.topoSort()){
                iscyc=-1;
            }
            else{
                iscyc=1;
                break;
            }
        }
        if(i==1){
            result=result+myGraph.myfunc+"\n";
            for(int k=0;k<variables.size();k++){
                dervs=dervs+"d"+myGraph.myfunc+"/"+"d"+variables[k]+" ";
            }
            dervs=dervs+"\n";
        }
        for(int x=0;x<myGraph.nodes.size();x++){
            if(myGraph.nodes[x]->name.compare(myGraph.myfunc)==0){
                stringstream ss;
                string out;
                ss<<myGraph.nodes[x]->myVal;
                ss>>out;
                result = result + out;
                result = result+"\n";
            }
        }
        for(int k=0;k<variables.size();k++){
            stringstream ss;
            string out;
            ss<<myGraph.nodes[myGraph.nodes.size()-1]->calculateDiff(variables[k]);
            ss>>out;
            dervs=dervs+out+" ";
        }
        dervs=dervs+"\n";
        for(int k=0;k<myGraph.nodes.size();k++){
            delete myGraph.nodes[k];
        }
    }
    if(iscyc==-1){
        ofstream rezs;
        rezs.open(argv[3]);
        rezs<<result;
        rezs.close();

        rezs.open(argv[4]);
        rezs<<dervs;
        rezs.close();
    }
    if(iscyc==1){
        cout<<"ERROR: COMPUTATION GRAPH HAS A CYCLE";
    }
}