#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include "Graph.h"
#include "Problem.h"
#include <ctime>
#include <queue>
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
    // "latin_square_10.col",
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
    "anna.col",
    "david.col",
    "games120.col",
    "homer.col",
    "huck.col",
    "jean.col",
    "miles1000.col",
    "miles1500.col",
    "miles250.col",
    "miles500.col",
    "miles750.col",
    "queen10_10.col",
    "queen11_11.col",
    "queen12_12.col",
    "queen13_13.col",
    "queen14_14.col",
    "queen15_15.col",
    "queen16_16.col",
    "queen5_5.col",
    "queen6_6.col",
    "queen7_7.col",
    "queen8_12.col",
    "queen8_8.col",
    "queen9_9.col",
};
const vector<string> MYC_FILES = {
    "myciel3.col",
    "myciel4.col",
    "myciel5.col",
    "myciel6.col",
    "myciel7.col",
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

void processaInstancia(vector<string> namesFiles, vector<DataGraph> &graphsData, double *porcentagem_processamento, double porcentagem_each_graph)
{

    cout << "Comecando" << endl;

    Problem problem;
    for (int i = 0; i < namesFiles.size(); i++)
    {
        cout << "Iniciando leitura do arquivo " << namesFiles[i] << endl;
        Graph graph;
        graph.read_from_file(namesFiles[i]);
        HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(graph);
        DataGraph dataGraph;
        dataGraph.nameGraph = namesFiles[i];
        dataGraph.data = response;
        graphsData.push_back(dataGraph);
        cout << "-- Grafo [" << namesFiles[i] << "]com nVertices: " << graph.get_num_vertices();
        cout << " - Cores: " << response.num_cores_encontradas
             << " - Conflitos: " << response.num_conflitos
             << " - Clique: " << response.size_clique
             << " - Delta: " << response.delta_clique
             << " - Arestas: " << response.num_arestas
             << " - DArestas: " << response.densidade_arestas << endl;
        cout << " - Cores genético: " << response.num_cores_genetico << " - Conflitos genético: " << response.num_conflitos_genetico << endl;
        cout << "Finalizando leitura do arquivo " << namesFiles[i] << endl;
        *porcentagem_processamento += porcentagem_each_graph;
        cout << "Porcentagem de processamento: " << *porcentagem_processamento << "%" << endl;
    }
};

int main()
{
    auto inicio = high_resolution_clock::now();
    cout << "Comecando" << endl;
    cout << "Iniciando leitura do arquivo" << endl;

    int total_graphs = 37;
    total_graphs += DSJ_FILES.size();
    total_graphs += REG_FILES.size();
    total_graphs += LAT_FILES.size();
    total_graphs += LEI_FILES.size();
    total_graphs += SCH_FILES.size();
    total_graphs += SGB_FILES.size();
    total_graphs += MYC_FILES.size();
    total_graphs += MIZ_FILES.size();
    total_graphs += HOS_FILES.size();
    total_graphs += CAR_FILES.size();
    total_graphs += CUL_FILES.size();
    total_graphs += R_FILES.size();

    double porcentagem_each_graph = 100.0 / total_graphs;
    double porcentagem_processamento = 0.0;

    ifstream file("grafo_aleatorio.txt");
    if (!file.is_open())
    {
        cout << "Error opening file" << endl;
        return 1;
    }

    cout << "Arquivo aberto" << endl;

    vector<DataGraph> graphsData;

    string line;

    Problem problem;
    int actual_graph = -1;
    Graph *graph;
    while (getline(file, line))
    {
        // cout << "lendo linha" << actual_graph + 1 << endl;
        stringstream linestream;
        linestream << line;
        char info[10];
        linestream.getline(info, 4, ':');
        if (!line.empty() && atoi(info) != actual_graph)
        {
            if (actual_graph != -1)
            {
                graph->generate_adj_matrix();
                cout << "-- Grafo [" << actual_graph << "]com nVertices: " << graph->get_num_vertices() << endl;
                HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(*graph);
                DataGraph dataGraph;
                dataGraph.nameGraph = "Graph " + to_string(actual_graph);
                dataGraph.data = response;
                graphsData.push_back(dataGraph);
                cout << " - Cores: " << response.num_cores_encontradas
                     << " - Clique: " << response.size_clique
                     << " - Delta: " << response.delta_clique
                     << " - Arestas: " << response.num_arestas
                     << " - DArestas: " << response.densidade_arestas << endl;
                cout << " - Cores genético: " << response.num_cores_genetico << " - Conflitos genético: " << response.num_conflitos_genetico << endl;
                porcentagem_processamento += porcentagem_each_graph;
                cout << "Porcentagem de processamento: " << porcentagem_processamento << "%" << endl;
            }
            actual_graph++;
            graph = new Graph();
        }

        graph->get_from_line(line);
    }
    file.close();

    processaInstancia(DSJ_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(REG_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(LAT_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(LEI_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(SCH_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(SGB_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(MYC_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(MIZ_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(HOS_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(CAR_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(CUL_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    processaInstancia(R_FILES, graphsData, &porcentagem_processamento, porcentagem_each_graph);
    auto fim = high_resolution_clock::now();

    // Calcula a duração
    auto duracao = duration_cast<milliseconds>(fim - inicio);

    cout << "Tempo de execução total: " << duracao.count() << " ms" << endl;
    fstream fileDataGraphs;
    fileDataGraphs.open("graphsData.csv", ios::out);
    // fileDataGraphs << "name_graph,num_vertices,solution_colors,n_conflitos,n_cores_bl_pmtv,n_conflitos_bl_pmtv,n_cores_bl_mmtv,n_conflitos_bl_mmtv, n_cores_bl_pmcmf,n_conflitos_bl_pmcmf,n_cores_bl_mmcmf,n_conflitos_bl_mmcmf,n_cores_tabu_troca_vertice,n_conflitos_tabu_troca_vertice,n_cores_tabu_menor_frequencia,n_conflitos_tabu_menor_frequencia,\n";
    fileDataGraphs << "name_graph,num_vertices,solution_colors,n_conflitos,n_cores_genetico,n_conflitos_genetico,\n";
    for (int i = 0; i < graphsData.size(); i++)
    {
        fileDataGraphs << graphsData[i].nameGraph                         // name_graph
                       << "," << graphsData[i].data.num_vertices          // num_vertices
                       << "," << graphsData[i].data.num_arestas           // num_arestas
                       << "," << graphsData[i].data.num_cores_encontradas // solution_colors
                       << "," << graphsData[i].data.num_conflitos         // n_conflitos
                       //    << "," << graphsData[i].data.num_cores_busca_local_primeira_melhora_troca_vertice           // n_cores_bl_pmtv
                       //    << "," << graphsData[i].data.num_conflitos_busca_local_primeira_melhora_troca_vertice       // n_conflitos_bl_pmtv
                       //    << "," << graphsData[i].data.num_cores_busca_local_melhor_melhora_troca_vertice             // n_cores_bl_mmtv
                       //    << "," << graphsData[i].data.num_conflitos_busca_local_melhor_melhora_troca_vertice         // n_conflitos_bl_mmtv
                       //    << "," << graphsData[i].data.num_cores_busca_local_primeira_melhora_cor_menos_frequente     // n_cores_bl_pmcmf
                       //    << "," << graphsData[i].data.num_conflitos_busca_local_primeira_melhora_cor_menos_frequente // n_conflitos_bl_pmcmf
                       //    << "," << graphsData[i].data.num_cores_busca_local_melhor_melhora_cor_menos_frequente       // n_cores_bl_mmcmf
                       //    << "," << graphsData[i].data.num_conflitos_busca_local_melhor_melhora_cor_menos_frequente   // n_conflitos_bl_mmcmf
                       //    << "," << graphsData[i].data.num_cores_tabu_troca_vertice                                   // n_cores_tabu_troca_vertice
                       //    << "," << graphsData[i].data.num_conflitos_tabu_troca_vertice                               // n_conflitos_tabu_troca_vertice
                       //    << "," << graphsData[i].data.num_cores_tabu_menor_frequencia                                //  n_cores_tabu_menor_frequencia
                       //    << "," << graphsData[i].data.num_conflitos_tabu_menor_frequencia                            // n_conflitos_tabu_menor_frequencia
                       << "," << graphsData[i].data.num_cores_genetico     // n_cores_genetico
                       << "," << graphsData[i].data.num_conflitos_genetico // n_conflitos_genetico
                       << ",\n";
    }
    fileDataGraphs.close();
    return 0;
}
