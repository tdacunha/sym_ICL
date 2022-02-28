package lib

import (
	"log"
	"runtime"
)

func MemoryUsage() {
	m := runtime.MemStats{ }
	log.Println("Memory Usage:")

	runtime.ReadMemStats(&m)
	log.Printf("    In use: %5.2f GB, Total use: %5.2f GB, Idle: %5.2g GB",
		float64(m.Alloc)/1e9, float64(m.TotalAlloc)/1e9,
		float64(m.HeapIdle)/1e9)
	log.Printf("    Next GC cycle: %5.2 GB, GC pause time: %.2f s",
		float64(m.NextGC)/1e9, float64(m.PauseTotalNs)*1e9)
}
