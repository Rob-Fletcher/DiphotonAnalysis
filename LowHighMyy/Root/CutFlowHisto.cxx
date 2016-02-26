#include <iostream>
#include <iomanip>

#include <LowHighMyy/CutFlowHisto.h>

ClassImp(CutFlowHisto)

using namespace std;


CutFlowHisto::CutFlowHisto(const char* name, const char* title, int n_cuts)
    : TH1F(name, title, n_cuts, 0, n_cuts) { };


void CutFlowHisto::add_cut(std::string name, bool * ptr_cut)
{
  m_pass_flags.push_back( std::make_pair( name, ptr_cut) );
  GetXaxis()->SetBinLabel(m_pass_flags.size(), name.c_str());
}

void CutFlowHisto::add_cuts(std::vector< std::pair<std::string, bool> > &pass_flags){
  m_pass_flags.reserve(100);
  unsigned int size = pass_flags.size();
  for( unsigned int i = 0; i < size; i++ ){
    add_cut(pass_flags[i].first, &pass_flags[i].second);
    cout<<"Setting cut "<<pass_flags[i].first<<" "<<pass_flags[i].second<<" pointer: "<<&pass_flags[i].second<<endl;
  }
}

bool CutFlowHisto::event_accepted(){
  adjust_cutflow();
  return *m_pass_flags.back().second;
}

bool CutFlowHisto::pass_cut(std::string pass_cut){
  adjust_cutflow();
  bool passed(false);
  for( auto pass_flag : m_pass_flags ) if(pass_flag.first == pass_cut) passed = pass_flag.second;
  return passed;
}


void CutFlowHisto::fill_cutflow(std::string pass_name, bool passed)
{
/*
  std::cout<<"**** cutflow status\n";
  for( auto pass_flag : m_pass_flags ){
    cout<<pass_flag.first<<" "<<*pass_flag.second<<" pointer: "<<pass_flag.second<<endl;
  }
*/
  for( auto pass_flag : m_pass_flags ){
    if(pass_flag.first == pass_name){
      //std::cout<<"Setting flag "<<pass_name<<" to "<<passed<<endl;
      *pass_flag.second = passed;
    }
  }
  if(m_pass_flags.back().first == pass_name) adjust_cutflow();
}

void CutFlowHisto::finalize_cutflow()
{
  adjust_cutflow();
  unsigned int size = m_pass_flags.size();
  for( unsigned int i = 0; i < size; i++ ){
    if(*m_pass_flags[i].second) Fill(i+0.5);
  }
}

void CutFlowHisto::adjust_cutflow()
{
  unsigned int size = m_pass_flags.size();
  for( unsigned int i = 0; i < size; i++ ){
    if(!*m_pass_flags[i].second) if( (i+1) < m_pass_flags.size() ) *m_pass_flags[i+1].second = false;
  }
}

void CutFlowHisto::reset(){
  for( auto pass_flag : m_pass_flags )
    *pass_flag.second = false;
  *m_pass_flags[0].second = true;
}

void CutFlowHisto::print_cutflow() const
{
    cout << "**** Printing cuflow\n";
    for (int i = 1; i != GetNbinsX() + 1; ++i) {
	    const int n = GetBinContent(i);
	    if (n == 0) break;
	    cout << "\n "
	         << setw(25) << GetXaxis()->GetBinLabel(i)
	         << '\t'
	         << GetBinContent(i)
	         << '\t';
	    streamsize ss = cout.precision();
	    if (i == 1) cout << setw(10) << "    ";
	    else {
	        cout << setw(10) << setprecision(3) << (GetBinContent(i) - GetBinContent(i - 1)) / (float)(GetBinContent(i - 1)) * 100. << " %";
	    }
	    cout << setw(10) << setprecision(3) << (GetBinContent(i)) / (float)(GetBinContent(1)) * 100. << " %";
	cout.precision(ss);
    }
    cout << endl;
}

CutFlowHisto::~CutFlowHisto(){
  m_pass_flags.clear();
}
