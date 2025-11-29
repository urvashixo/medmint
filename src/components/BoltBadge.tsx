import React from 'react'

export function BoltBadge() {
  return (
    <a
      href="https://bolt.new"
      target="_blank"
      rel="noopener noreferrer"

      style={{ zIndex: 9999 }}
    >
      <img
        src="/white_circle_360x360.png"
        alt="Powered by Bolt"
        className="w-[50px] h-[50px] md:w-[100px] md:h-[100px] rounded-full shadow-lg hover:shadow-xl transition-shadow duration-300"
      />
    </a>
  )
}