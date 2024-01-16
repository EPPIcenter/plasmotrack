//
// Created by Maxwell Murphy on 6/20/23.
//

#include "core/io/parse_json.h"
#include "core/samplers/meta/ReplicaExchange.h"
#include "core/utils/timers.h"
#include "impl/model/ModelEight/Model.h"
#include "impl/model/ModelEight/ModelLogger.h"
#include "impl/model/ModelEight/SampleScheduler.h"
#include "impl/model/ModelEight/State.h"
#include "impl/model/ModelEight/StateLogger.h"

#include "Application.h"
#include <string>
#include <filesystem>
#include <vector>

#include <thread>


#include "imgui.h"

using namespace transmission_nets::impl;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::samplers;
using namespace transmission_nets::core::utils;
using namespace transmission_nets::impl;


namespace transmission_nets::tools::gui::ModelEightGUI {

    namespace fs = std::filesystem;
    struct TxNetGUI {

        std::string input = "/data/mmurphy/analysis/astmh/superinf_simple/nodes.json";
        std::string output_dir = "/data/mmurphy/analysis/Zanzibar_tx_nets/output_gui/";
        std::string symptomatic_idp_path = "/data/mmurphy/analysis/Zanzibar_tx_nets/symptomatic_idp.csv";
        std::string asymptomatic_idp_path = "/data/mmurphy/analysis/Zanzibar_tx_nets/asymptomatic_idp.csv";

        long seed = timers::time().time_since_epoch().count();
        bool null_model = false;
        int num_chains = 1;
        int num_cores = 1;
        double gradient = 1.0;

        bool hotload = false;


        const fs::path nodesFile{input};
        const fs::path outputDir{output_dir};
        const std::vector<long double> symptomatic_idp = loadVectorFromFile<long double>(symptomatic_idp_path);
        const std::vector<long double> asymptomatic_idp = loadVectorFromFile<long double>(asymptomatic_idp_path);

        std::ifstream inputFile{nodesFile};
        const json j = loadJSON(inputFile);

        std::shared_ptr<boost::random::mt19937> r = std::make_shared<boost::random::mt19937>(seed);
        std::unique_ptr<ReplicaExchange<ModelEight::State, ModelEight::Model, ModelEight::SampleScheduler, ModelEight::ModelLogger, ModelEight::StateLogger>> repex = std::make_unique<ReplicaExchange<ModelEight::State, ModelEight::Model, ModelEight::SampleScheduler, ModelEight::ModelLogger, ModelEight::StateLogger>>(num_chains, 100000, gradient, r, outputDir, hotload, null_model, num_cores, j, symptomatic_idp, asymptomatic_idp);

        [[noreturn]] void compute() {
            while (true) {
                repex->sample();
            }
        }

    };

    TxNetGUI txnet_gui;

    void Compute() {
        txnet_gui.compute();
    }

    void RenderUI() {

            ImGui::Begin("Model Eight");

            static float values[1024] = {};
            static int values_offset = 0;
            static float refresh_time = 0.0f;
            static int total_samples = 0;

            refresh_time = ImGui::GetTime();

            while (refresh_time < ImGui::GetTime()) {
                auto val = txnet_gui.repex->hotValueThreadSafe();
                // sometimes the hot value is transiently -infinity, which is not a valid value for the plot
                if (val <= -std::numeric_limits<typeof(val)>::infinity()) {
                    val = values[values_offset-1];
                }
                values[values_offset] = val;
                values_offset = (values_offset+1) % IM_ARRAYSIZE(values);
                refresh_time += 1.0f/60.0f;
                total_samples++;
            }

            {
                float average = 0.0f;
                float sd = 0.0f;
                for (int n = 0; n < IM_ARRAYSIZE(values); n++)
                    average += values[n];
                average /= (float)IM_ARRAYSIZE(values);

                for (int n = 0; n < IM_ARRAYSIZE(values); n++)
                    sd += (values[n] - average) * (values[n] - average);
                sd /= (float)IM_ARRAYSIZE(values);
                sd = sqrt(sd);

                char overlay[128];
                sprintf(overlay, "average: %f, total %i", average, total_samples);
                ImGui::PlotLines("##plot", values, IM_ARRAYSIZE(values), values_offset, overlay, average - sd * 3, average + sd * 3, ImVec2(0,100.0f));
            }

            ImGui::End();



            ImGui::ShowDemoWindow();
            ImGui::ShowDebugLogWindow();

    };

}
