#ifndef PHYLIP32M10STATE_H
#define PHYLIP32M10STATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class phylip32_m10_state : public BaseFormatHandler {
public:

    explicit phylip32_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip32_m10";
        extension = "phy";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // PHYLIP32STATE_H
