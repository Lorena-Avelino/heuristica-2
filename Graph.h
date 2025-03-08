#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cstring>

using namespace std;

class CoupleVertice
{
public:
    unsigned vertice1;
    unsigned vertice2;
};

class Graph
{

protected:
    vector<vector<CoupleVertice>> adj_list;
    vector<vector<bool>> adj_matrix;

public:
    void get_from_line(string line);
    void generate_adj_matrix();
    vector<CoupleVertice> get_neighbors(int vertice);
    int get_num_vertices();
    bool ehVizinho(int vertice1, int vertice2);
    void print_graph();
    vector<vector<CoupleVertice>> get_adj_list() { return adj_list; };
    void read_from_file(string filename);
};

void Graph::generate_adj_matrix()
{
    adj_matrix.resize(adj_list.size());
    for (int i = 0; i < adj_list.size(); i++)
    {
        adj_matrix[i].resize(adj_list.size());
        for (int j = 0; j < adj_list.size(); j++)
        {
            adj_matrix[i][j] = false;
        }
    }

    for (int i = 0; i < adj_list.size(); i++)
    {
        for (int j = 0; j < adj_list[i].size(); j++)
        {
            adj_matrix[i][adj_list[i][j].vertice2] = true;
        }
    }
}

bool Graph::ehVizinho(int vertice1, int vertice2)
{
    return adj_matrix[vertice1][vertice2] && adj_matrix[vertice2][vertice1];
}

void Graph::get_from_line(string line)
{
    vector<CoupleVertice> neighbors;

    stringstream linestream;
    linestream << line;
    char info[10];
    linestream.getline(info, 5, ':');
    linestream.getline(info, 5, ':');
    int vertice1 = atoi(info);

    int vertice2;
    char comma;
    while (linestream >> vertice2)
    {
        CoupleVertice couple;
        couple.vertice1 = vertice1;
        couple.vertice2 = vertice2;
        neighbors.push_back(couple);
        linestream >> comma;
    }

    if (adj_list.size() < vertice1 + 1)
    {
        adj_list.resize(vertice1 + 1);
    }

    adj_list[vertice1] = neighbors;
};

vector<CoupleVertice> Graph::get_neighbors(int vertice)
{
    return adj_list[vertice];
};

void Graph::print_graph()
{
    for (int i = 0; i < adj_list.size(); i++)
    {

        cout << "MAIN Vertex " << i << endl;

        for (int j = 0; j < adj_list[i].size(); j++)
        {
            cout << "Vertex " << adj_list[i][j].vertice1 << ":";
            cout << " " << adj_list[i][j].vertice2;
            cout << endl;
        }
    }
};

int Graph::get_num_vertices()
{
    return adj_list.size();
}

void Graph::read_from_file(string filename)
{
    ifstream file;
    file.open(filename);
    if (!file.good())
    {
        cout << "Erro ao abrir arquivo" << endl;
        exit(1);
    }

    string line;

    while (getline(file, line))
    {
        if (line[0] != 'c')
        {
            if (line[0] == 'p')
            {
                stringstream linestream;
                linestream << line;
                char info[10];
                linestream.getline(info, 5, ' ');
                linestream.getline(info, 5, ' ');
                linestream.getline(info, 5, ' ');
                int quantidade = atoi(info);
                adj_list.resize(quantidade);
            }
            if (line[0] == 'e')
            {
                stringstream linestream;
                linestream << line;
                char info[30];
                linestream.getline(info, 5, ' ');
                linestream.getline(info, 5, ' ');
                int vertice1 = atoi(info);
                linestream.getline(info, 5, ' ');
                int vertice2 = atoi(info);
                CoupleVertice couple;
                couple.vertice1 = vertice1 - 1;
                couple.vertice2 = vertice2 - 1;
                adj_list[couple.vertice1].push_back(couple);
            }
        }
    }
    file.close();
    generate_adj_matrix();
}
#endif