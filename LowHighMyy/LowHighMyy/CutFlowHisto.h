#ifndef MyAnalysis_CutFlowHisto_H
#define MyAnalysis_CutFlowHisto_H

#include <string>
#include <vector>

#include <TH1F.h>

// TODO: must be mergeable

class CutFlowHisto : public TH1F
{
public:
    CutFlowHisto() { };
    CutFlowHisto(const char* name, const char* title, int n_cuts = 20);
    ~CutFlowHisto();

    void add_cuts(std::vector< std::pair<std::string, bool> > &pass_flags);
    void fill_cutflow(std::string pass_name, bool pass_flag);
    void finalize_cutflow();
    void print_cutflow() const;
    void reset();
    bool event_accepted();
    bool pass_cut(std::string pass_cut);
    ClassDef(CutFlowHisto, 1);

private:
    std::vector< std::pair<std::string, bool*> > m_pass_flags; //!

    void adjust_cutflow();
    void add_cut(std::string name, bool * ptr_cut);
};

#endif
