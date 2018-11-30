#include <genmap-impl.h>

#include <stdio.h>
//
// GenmapHeader: Create, Destroy
//
//int GenmapCreateHeader(GenmapHandle h, GenmapHeader *header) {
//  return h->CreateHeader(header);
//}
//
//int GenmapDestroyHeader(GenmapHandle h, GenmapHeader header) {
//  return h->DestroyHeader(header);
//}
////
//// GenmapElements: Create, Destroy
////
//int GenmapCreateElements(GenmapHandle h, GenmapElements *e) {
//  return h->CreateElements(e);
//}
//
//int GenmapDestroyElements(GenmapHandle h, GenmapElements e) {
//  return h->DestroyElements(e);
//}
//
//GenmapElements GenmapGetElements(GenmapHandle h) {
//  return h->GetElements(h);
//}
//
// Do File I/O in parallel
//
int GenmapRead(GenmapHandle h, void* data) {
  return h->Read(h, data);
}

int GenmapWrite(GenmapHandle h, char *fileNameBase) {
  return h->Write(h, fileNameBase);
}
