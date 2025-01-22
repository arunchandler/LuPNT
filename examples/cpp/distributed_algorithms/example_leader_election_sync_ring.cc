#include <lupnt/core/scheduler.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

using namespace std;
using namespace lupnt;

// Transmission class
struct Transmission {
  int id_;
};

class MyTransonder;

// Channel class
class Channel {
private:
  int id_;
  Real delay_ = 0.5;
  vector<shared_ptr<MyTransonder>> transponders_;

public:
  void Send(Real t, const Transmission &transm, const MyTransonder &sender);
  void Add(shared_ptr<MyTransonder> MyTransonder) { transponders_.push_back(MyTransonder); }
  vector<shared_ptr<MyTransonder>> GetMyTransonders() { return transponders_; }
};

// MyTransonder class
class MyTransonder {
private:
  shared_ptr<Channel> channel_;
  function<void(Real, const Transmission &)> receive_callback_;

public:
  MyTransonder(shared_ptr<Channel> channel) : channel_(channel) {}
  void SetReceiveCallback(function<void(Real, const Transmission &)> callback) {
    receive_callback_ = callback;
  }

  void Send(Real t, const Transmission &transm) { channel_->Send(t, transm, *this); }

  void Receive(Real t, const Transmission &transm) { receive_callback_(t, transm); }
};

void Channel::Send(Real t, const Transmission &transm, const MyTransonder &sender) {
  for (const auto &MyTransonder : transponders_) {
    if (MyTransonder.get() != &sender) {
      Scheduler::Schedule(t + delay_,
                          [MyTransonder, transm](Real t) { MyTransonder->Receive(t, transm); });
    }
  }
}

// Agent class
class RingAgent {
private:
  static int id_counter_;
  const int id_;
  vector<shared_ptr<MyTransonder>> transponders_;

public:
  RingAgent() : id_(id_counter_++), transponders_(vector<shared_ptr<MyTransonder>>()) {}
  int GetId() { return id_; }

  void Add(shared_ptr<MyTransonder> MyTransonder) { transponders_.push_back(MyTransonder); }
  vector<shared_ptr<MyTransonder>> GetMyTransonders() { return transponders_; }
};

int RingAgent::id_counter_ = 0;

// Application class
class LeaderElectionSyncRingApp : public Application {
private:
  shared_ptr<RingAgent> agent_;
  int id_;
  int id_received_ = -1;
  bool is_leader_ = false;

public:
  LeaderElectionSyncRingApp(shared_ptr<RingAgent> agent) : agent_(agent), id_(agent->GetId()) {}
  Real GetFrequency() override { return 1.0; }
  void Step(Real t) override {
    if (id_received_ == -1) {
      cout << "[Agent " << id_ << "] Sending " << id_ << " at t = " << t << endl;
      for (const auto &MyTransonder : agent_->GetMyTransonders()) {
        // transmittion instance with _id
        MyTransonder->Send(t, Transmission{id_});
      }
    } else if (id_received_ > id_) {
      cout << "[Agent " << id_ << "] Sending " << id_received_ << " at t = " << t << endl;
      for (const auto &MyTransonder : agent_->GetMyTransonders()) {
        MyTransonder->Send(t, Transmission{id_received_});
      }
    } else if (id_received_ < id_) {
      cout << "[Agent " << id_ << "] Doing nothing"
           << " at t = " << t << endl;
    } else {
      is_leader_ = true;
      cout << "[Agent " << id_ << "] Leader"
           << " at t = " << t << endl;
    }
  }
  void TransmissionReceived(Real t, const Transmission &transm) {
    cout << "[Agent " << id_ << "] Received " << transm.id_ << " at t = " << t << endl;
    id_received_ = max(id_received_, transm.id_);
  }
};

int main() {
  int n = 3;
  vector<shared_ptr<RingAgent>> agents(n);
  vector<shared_ptr<LeaderElectionSyncRingApp>> apps(n);
  vector<shared_ptr<Channel>> channels(n);

  // Create agents, applications, and channels
  for (int i = 0; i < n; i++) {
    agents[i] = make_shared<RingAgent>();
    apps[i] = make_shared<LeaderElectionSyncRingApp>(agents[i]);
    channels[i] = make_shared<Channel>();
  }

  // Create MyTransonders and add them to agents and channels
  for (int i = 0; i < n; i++) {
    auto ch1 = channels[i];
    auto ch2 = channels[(i + 1) % n];
    auto transc1 = make_shared<MyTransonder>(ch1);
    auto transc2 = make_shared<MyTransonder>(ch2);
    transc1->SetReceiveCallback(bind(&LeaderElectionSyncRingApp::TransmissionReceived,
                                     apps[i].get(), placeholders::_1, placeholders::_2));
    transc2->SetReceiveCallback(bind(&LeaderElectionSyncRingApp::TransmissionReceived,
                                     apps[i].get(), placeholders::_1, placeholders::_2));
    agents[i]->Add(transc1);
    agents[i]->Add(transc2);
    ch1->Add(transc1);
    ch2->Add(transc2);
  }

  // Schedule applications
  Real t_start = 0.0;
  Real freq = 1.0;
  for (int i = 0; i < n; i++) Scheduler::ScheduleApplication(*apps[i], t_start);

  // Run simulation
  Scheduler::RunSimulation(4.0);
  return 0;
}
