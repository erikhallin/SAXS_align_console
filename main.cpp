#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <cmath>

const float _pi=3.14159265359;

using namespace std;

struct st_pos
{
    st_pos()
    {
        x=y=0.0;
    }
    st_pos(float _x,float _y,float _z)
    {
        x=_x;
        y=_y;
        z=_z;
    }
    /*st_pos(st_pos& _pos)
    {
        x=_pos.x;
        y=_pos.y;
    }*/

    float x,y,z;

    st_pos operator-(st_pos _pos)
    {
        return st_pos(x-_pos.x,y-_pos.y,z-_pos.z);
    }
    st_pos operator-=(st_pos _pos)
    {
        x-=_pos.x;
        y-=_pos.y;
        z-=_pos.z;

        return *this;
    }
    st_pos operator+(st_pos _pos)
    {
        return st_pos(x+_pos.x,y+_pos.y,z+_pos.z);
    }
    st_pos operator+=(st_pos _pos)
    {
        x+=_pos.x;
        y+=_pos.y;
        z+=_pos.z;

        return *this;
    }
    st_pos operator*(float scale)
    {
        return st_pos(x*scale,y*scale,z*scale);
    }
    st_pos operator*=(float scale)
    {
        x*=scale;
        y*=scale;
        z*=scale;

        return *this;
    }
    st_pos operator/(float scale)
    {
        return st_pos(x/scale,y/scale,z/scale);
    }
    st_pos operator/=(float scale)
    {
        x/=scale;
        y/=scale;
        z/=scale;

        return *this;
    }
    bool operator==(st_pos _pos)
    {
        return (x==_pos.x && y==_pos.y && z==_pos.z);
    }
    float distance3(st_pos _pos)
    {
        return ( (x-_pos.x)*(x-_pos.x)+(y-_pos.y)*(y-_pos.y)+(z-_pos.z)*(z-_pos.z) );
    }
    float distance(st_pos _pos)
    {
        return sqrt(distance3(_pos));
    }
    float length(void)
    {
        return sqrt( x*x+y*y+z*z );
    }
    float normalize()
    {
        float length=sqrt( x*x+y*y+z*z );
        x/=length;
        y/=length;
        z/=length;

        return length;
    }

};

string float_to_pdb_string(float value);

int main()
{
    cout<<"Bead model alignment\n\n";

    float radius_cutoff=3.0;
    vector<st_pos> vec_model1;
    vector<st_pos> vec_model2;
    st_pos center_model1;
    st_pos center_model2;

    //load models
    cout<<"Name of File 1: ";
    string input,input_file_name1,input_file_name2;
    while(true)
    {
        getline(cin,input);

        //test if file exists
        ifstream file(input.c_str());
        if(file==0)
        {
            cout<<"\nERROR: Could not find that file\n";
            file.close();
        }
        else
        {
            file.close();
            break;
        }
    }
    input_file_name1=input;
    cout<<"Name of File 2: ";
    while(true)
    {
        getline(cin,input);

        //test if file exists
        ifstream file(input.c_str());
        if(file==0)
        {
            cout<<"\nERROR: Could not find that file\n";
            file.close();
        }
        else
        {
            file.close();
            break;
        }
    }
    input_file_name2=input;

    //get bead positions
    ifstream input_file1(input_file_name1.c_str());
    string word,line;
    float pos_avg[3]={0,0,0};
    int atom_counter=0;
    while(getline(input_file1,line))
    {
        //copy if not starts with ATOM
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            atom_counter++;

            float pos[3]={0,0,0};
            //read and scale position
            pos[0]=atof(string(line,30,8).c_str());
            pos[1]=atof(string(line,38,8).c_str());
            pos[2]=atof(string(line,46,8).c_str());

            pos_avg[0]+=pos[0];
            pos_avg[1]+=pos[1];
            pos_avg[2]+=pos[2];

            vec_model1.push_back(st_pos(pos[0],pos[1],pos[2]));
        }
    }
    input_file1.close();
    pos_avg[0]/=(float)atom_counter;
    pos_avg[1]/=(float)atom_counter;
    pos_avg[2]/=(float)atom_counter;
    center_model1=st_pos(pos_avg[0],pos_avg[1],pos_avg[2]);
    //move to center
    for(unsigned int beadi=0;beadi<vec_model1.size();beadi++)
    {
        vec_model1[beadi]-=center_model1;
    }

    ifstream input_file2(input_file_name2.c_str());
    pos_avg[0]=pos_avg[1]=pos_avg[2]=0;
    atom_counter=0;
    while(getline(input_file2,line))
    {
        //copy if not starts with ATOM
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            atom_counter++;

            float pos[3]={0,0,0};
            //read and scale position
            pos[0]=atof(string(line,30,8).c_str());
            pos[1]=atof(string(line,38,8).c_str());
            pos[2]=atof(string(line,46,8).c_str());

            pos_avg[0]+=pos[0];
            pos_avg[1]+=pos[1];
            pos_avg[2]+=pos[2];

            vec_model2.push_back(st_pos(pos[0],pos[1],pos[2]));
        }
    }
    input_file2.close();
    pos_avg[0]/=(float)atom_counter;
    pos_avg[1]/=(float)atom_counter;
    pos_avg[2]/=(float)atom_counter;
    center_model2=st_pos(pos_avg[0],pos_avg[1],pos_avg[2]);
    //move to center
    for(unsigned int beadi=0;beadi<vec_model2.size();beadi++)
    {
        vec_model2[beadi]-=center_model2;
    }

    //overlap test
    int beads_inside=0;
    int beads_outside=0;
    float min_dist3=radius_cutoff*radius_cutoff*radius_cutoff;
    for(unsigned int beadi1=0;beadi1<vec_model1.size();beadi1++)
    {
        bool bead_inside=false;
        for(unsigned int beadi2=0;beadi2<vec_model2.size();beadi2++)
        {
            //test distance
            if( vec_model1[beadi1].distance3(vec_model2[beadi2])<min_dist3 )
            {
                bead_inside=true;
                //is within other model
                beads_inside++;
                break;
            }

        }

        //bead was not inside other model
        if(!bead_inside) beads_outside++;
    }

    //cout<<"1 center: "<<center_model2.x<<", "<<center_model2.y<<", "<<center_model2.z<<"\n";
    cout<<(int)vec_model1.size()<<"\tbeads in model 1\n";
    cout<<(int)vec_model2.size()<<"\tbeads in model 2\n";
    cout<<beads_inside<<"\tBeads inside\n";
    cout<<beads_outside<<"\tBeads outside\n";


    //position optimisation, move model 1 to fit model 2
    float high_score=0;
    int cycle_counter=0;
    while(true)
    {
        cycle_counter++;

        //translate
        float trans_sens=4.0;
        int rand_acc=10000;
        float x_trans=trans_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        float y_trans=trans_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        float z_trans=trans_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        for(unsigned int beadi=0;beadi<vec_model1.size();beadi++)
        {
            vec_model2[beadi].x+=x_trans;
            vec_model2[beadi].y+=y_trans;
            vec_model2[beadi].z+=z_trans;
        }

        //rotate
        float rot_sens=_pi*0.5;
        float x_rot=rot_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        float y_rot=rot_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        float z_rot=rot_sens*((rand()%rand_acc)/(float)rand_acc*2.0-1.0);
        float sin_x = sin(x_rot);
        float cos_x = cos(x_rot);
        float sin_y = sin(y_rot);
        float cos_y = cos(y_rot);
        float sin_z = sin(z_rot);
        float cos_z = cos(z_rot);
        float x_pos,y_pos,z_pos;
        for(unsigned int beadi=0;beadi<vec_model1.size();beadi++)
        {
            //z
            x_pos=vec_model1[beadi].x;
            y_pos=vec_model1[beadi].y;
            vec_model1[beadi].x=x_pos*cos_z-y_pos*sin_z;
            vec_model1[beadi].y=y_pos*cos_z+x_pos*sin_z;
            //y
            x_pos=vec_model1[beadi].x;
            z_pos=vec_model1[beadi].z;
            vec_model1[beadi].x=x_pos*cos_y-z_pos*sin_y;
            vec_model1[beadi].z=z_pos*cos_y+x_pos*sin_y;
            //x
            y_pos=vec_model1[beadi].y;
            z_pos=vec_model1[beadi].z;
            vec_model1[beadi].y=y_pos*cos_x-z_pos*sin_x;
            vec_model1[beadi].z=z_pos*cos_x+y_pos*sin_x;
        }

        //score
        beads_inside=0;
        beads_outside=0;
        for(unsigned int beadi1=0;beadi1<vec_model1.size();beadi1++)
        {
            bool bead_inside=false;
            for(unsigned int beadi2=0;beadi2<vec_model2.size();beadi2++)
            {
                //test distance
                if( vec_model1[beadi1].distance3(vec_model2[beadi2])<min_dist3 )
                {
                    bead_inside=true;
                    //is within other model
                    beads_inside++;
                    break;
                }

            }

            //bead was not inside other model
            if(!bead_inside) beads_outside++;
        }

        //break test
        float score=(float)beads_inside/(int)vec_model1.size();
        if(score>high_score)
        {
            high_score=score;
            cout<<high_score<<endl;
        }
        if(cycle_counter>1000)
        {
            cycle_counter=0;
            //ask if continue
            cout<<"Type [stop] to stop\n";
            string msg;
            getline(cin,msg);
            if(msg.length()>0) break;
        }
        break;
    }

    //make new model file 1
    ofstream new_file1("new_file1.pdb");
    if(new_file1==0)
    {
        cout<<"ERROR: Could not create output file\n";
        return 1;
    }
    input_file1.open(input_file_name1.c_str());
    atom_counter=0;
    while(getline(input_file1,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            cout<<"start ";
            if(atom_counter>=(int)vec_model1.size())
            {
                cout<<"ERROR: Bad atom number\n";
                break;
            }
            cout<<"next ";

            //replace with new values
            string x_string=float_to_pdb_string( vec_model1[atom_counter].x );
            string y_string=float_to_pdb_string( vec_model1[atom_counter].y );
            string z_string=float_to_pdb_string( vec_model1[atom_counter].z );
            //cout<<x_string<<y_string<<z_string<<endl;
            cout<<"sec ";
            int line_pos_x=30;
            int line_pos_y=38;
            int line_pos_z=46;
            for(int i=0;i<8;i++)
            {
                line[line_pos_x+i]=x_string[i];
                line[line_pos_y+i]=y_string[i];
                line[line_pos_z+i]=z_string[i];
            }

            atom_counter++;
        }
        //cout<<line<<endl;

        //line to new file
        //new_file1<<line<<endl;
    }
    new_file1.close();

    //make new model file 2
    ofstream new_file2("new_file2.pdb");
    if(new_file2==0)
    {
        cout<<"ERROR: Could not create output file\n";
        return 1;
    }
    input_file2.open(input_file_name2.c_str());
    atom_counter=0;
    while(getline(input_file2,line))
    {
        if(line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
        {
            if(atom_counter>=(int)vec_model2.size())
            {
                cout<<"ERROR: Bad atom number\n";
                break;
            }

            //replace with new values
            string x_string=float_to_pdb_string( vec_model2[atom_counter].x );
            string y_string=float_to_pdb_string( vec_model2[atom_counter].y );
            string z_string=float_to_pdb_string( vec_model2[atom_counter].z );
            int line_pos_x=30;
            int line_pos_y=38;
            int line_pos_z=46;
            for(int i=0;i<8;i++)
            {
                line[line_pos_x+i]=x_string[i];
                line[line_pos_y+i]=y_string[i];
                line[line_pos_z+i]=z_string[i];
            }

            atom_counter++;
        }

        //line to new file
        new_file2<<line<<endl;
    }
    new_file2.close();

    cout<<"\nComplete\n\n";

    return 0;
}

string float_to_pdb_string(float value)
{
    //return string("   1.000");

    if(value>=10000.000) return string("  ERROR ");
    if(value<=-1000.000) return string("  ERROR ");

    //rounding
    if(value>0) value+=0.0005;
    if(value<0) value-=0.0005;

    //value convert to strings and reformat to length 8
    stringstream ss;
    ss<<(float)value;
    string string_value(ss.str());
    bool is_float=false;
    for(int i=0;i<string_value.length();i++)
    {
        if(string_value[i]=='.')
        {
            is_float=true;
            //cut if more than 3 decimals
            if((int)string_value.length()-i>4)
            {
                cout<<"cut1 "<<string_value<<" "<<(int)string_value.length()<<" "<<i;

                string cut_string=string(string_value,0,i+4);//crash here XXXXXXXXXXXXXXXXXXXXXXXXX
                string_value=cut_string;
                cout<<"cut2 "<<string_value<<endl;
            }
            break;
        }
    }
    if(!is_float) string_value.append(".000");
    string ret_string(8,' ');
    int counter=0;
    for(int i=7;i>=0;i--)
    {
        if((int)string_value.size()-1-counter<0)
        {
            //space fill
            ret_string[i]=' ';
        }
        else ret_string[i]=string_value[ (int)string_value.size()-1-counter ];
        counter++;
    }

    return ret_string;
}
