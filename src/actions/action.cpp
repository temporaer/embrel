#include <stdexcept> 
#include "action.hpp"
#include <factory/factory.h>
#include <nana.h>

using namespace std;

void Action::operator()()
{
	throw logic_error("Called Action() w/o subclassing");
}
void Action::configure()
{
}
Action::~Action()
{
}

namespace{ registerInFactory<Action, Action> registerBase("Action"); }
