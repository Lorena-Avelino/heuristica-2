#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include "Graph.h"
#include "Problem.h"
#include <ctime>
#include <set>
#include <algorithm>

using namespace std;
const vector<string> DSJ_FILES = {
    "DSJC1000.1.col",
    "DSJC1000.5.col",
    "DSJC1000.9.col",
    "DSJC125.1.col",
    "DSJC125.5.col",
    "DSJC125.9.col",
    "DSJC250.1.col",
    "DSJC250.5.col",
    "DSJC250.9.col",
    "DSJC500.1.col",
    "DSJC500.5.col",
    "DSJC500.9.col",
    "DSJR500.1.col",
    "DSJR500.1c.col",
    "DSJR500.5.col",
};
const vector<string> REG_FILES = {
    "fpsol2.i.1.col",
    "fpsol2.i.2.col",
    "fpsol2.i.3.col",
    "inithx.i.1.col",
    "inithx.i.2.col",
    "inithx.i.3.col",
    "mulsol.i.1.col",
    "mulsol.i.2.col",
    "mulsol.i.3.col",
    "mulsol.i.4.col",
    "mulsol.i.5.col",
    "zeroin.i.1.col",
    "zeroin.i.2.col",
    "zeroin.i.3.col",
};
const vector<string> LAT_FILES = {
    "latin_square_10.col",
};
const vector<string> LEI_FILES = {
    "le450_15a.col",
    "le450_15b.col",
    "le450_15b.col.txt",
    "le450_15c.col",
    "le450_15d.col",
    "le450_25a.col",
    "le450_25b.col",
    "le450_25c.col",
    "le450_25d.col",
    "le450_5a.col",
    "le450_5b.col",
    "le450_5c.col",
    "le450_5d.col",
};
const vector<string> SCH_FILES = {
    "school1.col",
    "school1_nsh.col",
};
const vector<string> SGB_FILES = {
    // "anna.col",
    // "david.col",
    // "games120.col",
    // "homer.col",
    // "huck.col",
    // "jean.col",
    // "miles1000.col",
    // "miles1500.col",
    // "miles250.col",
    // "miles500.col",
    // "miles750.col",
    // "queen10_10.col",
    // "queen11_11.col",
    // "queen12_12.col",
    // "queen13_13.col",
    // "queen14_14.col",
    // "queen15_15.col",
    // "queen16_16.col",
    "queen5_5.col",
    // "queen6_6.col",
    // "queen7_7.col",
    // "queen8_12.col",
    // "queen8_8.col",
    // "queen9_9.col",
};
const vector<string> MYC_FILES = {
    "myciel3.col",
    "myciel4.col",
    // "myciel5.col",
    // "myciel6.col",
    // "myciel7.col",
};
const vector<string> MIZ_FILES = {
    "mug100_1.col",
    "mug100_25.col",
    "mug88_1.col",
    "mug88_25.col",
};
const vector<string> HOS_FILES = {};
const vector<string> CAR_FILES = {
    "1-FullIns_3.col",
    "1-FullIns_4.col",
    "1-FullIns_5.col",
    "1-Insertions_4.col",
    "1-Insertions_5.col",
    "1-Insertions_6.col",
    "2-FullIns_3.col",
    "2-FullIns_4.col",
    "2-FullIns_5.col",
    "2-Insertions_3.col",
    "2-Insertions_4.col",
    "2-Insertions_5.col",
    "3-FullIns_3.col",
    "3-FullIns_4.col",
    "3-FullIns_5.col",
    "3-Insertions_3.col",
    "3-Insertions_4.col",
    "3-Insertions_5.col",
    "4-FullIns_3.col",
    "4-FullIns_4.col",
    "4-FullIns_5.col",
    "4-Insertions_3.col",
    "4-Insertions_4.col",
    "5-FullIns_3.col",
    "5-FullIns_4.col",
};
const vector<string> CUL_FILES = {
    "flat300_28_0.col",
    "flat1000_76_0.col",
    "flat1000_60_0.col",
    "flat1000_50_0.col",
};
const vector<string> R_FILES = {
    "r250.1c.col",
    "r250.5.col",
    "r1000.1c.col",
    "r1000.5.col",
};

struct DataGraph
{
    string nameGraph;
    HeuristicaConstrutivaResponse data;
};

void processaInstancia(vector<string> namesFiles, vector<DataGraph> &graphsData)
{
    Problem problem;
    for (int i = 0; i < namesFiles.size(); i++)
    {
        Graph graph;
        graph.read_from_file(namesFiles[i]);
        HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(graph);
        DataGraph dataGraph;
        dataGraph.nameGraph = namesFiles[i];
        dataGraph.data = response;
        graphsData.push_back(dataGraph);
        cout << "-- Grafo [" << namesFiles[i] << "]com nVertices: " << graph.get_num_vertices();
        cout << " - Cores: " << response.num_cores_encontradas
             << " - Clique: " << response.size_clique
             << " - Delta: " << response.delta_clique
             << " - Arestas: " << response.num_arestas
             << " - DArestas: " << response.densidade_arestas << endl
             << "------------------------------------------------------------------------------"
             << endl;
    }
};

vector<int> substitui_cor_menos_frequente(vector<int> &cores, Graph &graph)
{
    int num_vertices = graph.get_num_vertices();
    unordered_map<int, int> frequencias;

    for (int cor : cores)
    {
        frequencias[cor]++;
    }

    int cor_menos_frequente = -1;
    int menor_frequencia = num_vertices + 1;
    for (auto &par : frequencias)
    {
        if (par.second < menor_frequencia)
        {
            menor_frequencia = par.second;
            cor_menos_frequente = par.first;
        }
    }

    for (int i = 0; i < num_vertices; ++i)
    {
        if (cores[i] == cor_menos_frequente)
        {
            set<int> cores_vizinhas;
            for (CoupleVertice vizinho : graph.get_neighbors(i))
            {
                cores_vizinhas.insert(cores[vizinho.vertice2]);
            }

            for (int nova_cor = 1;; ++nova_cor)
            {
                if (cores_vizinhas.find(nova_cor) == cores_vizinhas.end())
                {
                    cores[i] = nova_cor;
                    break;
                }
            }
        }
    }

    return cores;
}

vector<int> busca_local(Graph &graph, vector<int> &solucao_inicial, bool melhorMelhora, vector<int> (*vizinhanca)(vector<int> &, Graph &))
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_solucao = solucao;
    int melhor_custo = *max_element(solucao.begin(), solucao.end());

    while (true)
    {
        vector<int> vizinho = vizinhanca(solucao, graph);
        int custo_vizinho = *max_element(vizinho.begin(), vizinho.end());

        if (custo_vizinho < melhor_custo)
        {
            melhor_solucao = vizinho;
            melhor_custo = custo_vizinho;

            if (!melhorMelhora)
            {
                break;
            }
        }
        else if (!melhorMelhora)
        {
            break;
        }
    }

    return melhor_solucao;
}

int main()
{

    ifstream file("grafo_aleatorio.txt");
    if (!file.is_open())
    {
        cout << "Error opening file" << endl;
        return 1;
    }

    vector<DataGraph> graphsData;

    string line;

    Problem problem;
    int actual_graph = -1;
    Graph *graph;
    // while (getline(file, line))
    // {
    //     stringstream linestream;
    //     linestream << line;
    //     char info[10];
    //     linestream.getline(info, 4, ':');
    //     if (!line.empty() && atoi(info) != actual_graph)
    //     {
    //         if (actual_graph != -1)
    //         {
    //             graph->generate_adj_matrix();
    //             cout << "-- Grafo [" << actual_graph << "]com nVertices: " << graph->get_num_vertices();
    //             HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(*graph);
    //             DataGraph dataGraph;
    //             dataGraph.nameGraph = "Graph " + to_string(actual_graph);
    //             dataGraph.data = response;
    //             graphsData.push_back(dataGraph);
    //             cout << " - Cores: " << response.num_cores_encontradas
    //                  << " - Clique: " << response.size_clique
    //                  << " - Delta: " << response.delta_clique
    //                  << " - Arestas: " << response.num_arestas
    //                  << " - DArestas: " << response.densidade_arestas << endl;
    //         }
    //         actual_graph++;
    //         graph = new Graph();
    //     }

    //     graph->get_from_line(line);
    // }
    file.close();

    processaInstancia(DSJ_FILES, graphsData);
    processaInstancia(REG_FILES, graphsData);
    processaInstancia(LAT_FILES, graphsData);
    processaInstancia(LEI_FILES, graphsData);
    processaInstancia(SCH_FILES, graphsData);
    processaInstancia(SGB_FILES, graphsData);
    processaInstancia(MYC_FILES, graphsData);
    processaInstancia(MIZ_FILES, graphsData);
    processaInstancia(HOS_FILES, graphsData);
    processaInstancia(CAR_FILES, graphsData);
    processaInstancia(CUL_FILES, graphsData);
    processaInstancia(R_FILES, graphsData);

    fstream fileDataGraphs;
    fileDataGraphs.open("graphsData.csv", ios::out);
    fileDataGraphs << "name_graph,num_vertices,num_arestas,clique_size,solution_colors,\n";
    for (int i = 0; i < graphsData.size(); i++)
    {
        fileDataGraphs << graphsData[i].nameGraph << "," << graphsData[i].data.num_vertices << "," << graphsData[i].data.num_arestas << "," << graphsData[i].data.size_clique << "," << graphsData[i].data.num_cores_encontradas << ",\n";
    }
    fileDataGraphs.close();
    return 0;
}
