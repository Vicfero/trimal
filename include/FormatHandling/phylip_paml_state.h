#ifndef PHYLIPPAMLSTATE_H
#define PHYLIPPAMLSTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class phylip_paml_state : public BaseFormatHandler {
public:

    explicit phylip_paml_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip_paml";
        extension = "phy";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // PHYLIPPAMLSTATE_H
