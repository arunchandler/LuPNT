#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  GroundStationDataMap gs_data = LoadGroundStationData();

  // print the ground station latitude, longitude, and altitude
  std::cout << "GS Name       Latitude    Longitude    Altitude" << std::endl;
  for (auto& [gs_id, gs] : gs_data) {
    std::cout << gs.name << "  " << gs.latitude << "  " << gs.longitude << "  " << gs.altitude_m
              << std::endl;
  }

  std::cout << " " << std::endl;
  GroundStation gs = GroundStation(gs_data[0]);
  std::cout << "GS Name  : " << gs.GetName() << std::endl;
  std::cout << "Lat, Lon : " << gs.GetLatitudeDouble() * DEG << " , "
            << gs.GetLongitudeDouble() * DEG << std::endl;
  std::cout << "PosVel   : " << gs.GetRvState()->GetVec().transpose() << std::endl;
}
